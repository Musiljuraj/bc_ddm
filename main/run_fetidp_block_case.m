function out = run_fetidp_block_case(cfg)
%RUN_FETIDP_BLOCK_CASE  End-to-end FETI-DP block run with diagnostics.
%
% Chain:
%   build_problem_data -> setup_fetidp -> solve_fetidp -> reconstruct_fetidp_solution
%   + compare vs direct global reference solve (Kff\Ff)
%
% Required fields:
%   cfg.n, cfg.nSubX, cfg.nSubY
%
% Optional fields:
%   cfg.f_handle    : RHS f(x,y) (default = @(x,y) 1)
%   cfg.tol         : PCG tolerance (default 1e-10)
%   cfg.maxit       : PCG max iters (default 500)
%   cfg.setup_opts  : opts passed to setup_fetidp (default struct('tol_rank',1e-12))
%   cfg.seed        : deterministic RNG seed (default 0)
%
% Output struct out contains data, lambda, stats, reconstruction and diagnostics.

  if nargin ~= 1 || ~isstruct(cfg)
    error('run_fetidp_block_case: cfg must be a struct.');
  end
  req = {'n','nSubX','nSubY'};
  for k = 1:numel(req)
    if ~isfield(cfg, req{k})
      error('run_fetidp_block_case: cfg.%s is required.', req{k});
    end
  end

  n     = cfg.n;
  nSubX = cfg.nSubX;
  nSubY = cfg.nSubY;

  f_handle = @(x,y) 1.0;
  if isfield(cfg,'f_handle') && ~isempty(cfg.f_handle)
    f_handle = cfg.f_handle;
  end
  if ~isa(f_handle,'function_handle')
    error('run_fetidp_block_case: cfg.f_handle must be a function handle.');
  end

  tol = 1e-10;
  if isfield(cfg,'tol') && ~isempty(cfg.tol)
    tol = cfg.tol;
  end
  maxit = 500;
  if isfield(cfg,'maxit') && ~isempty(cfg.maxit)
    maxit = cfg.maxit;
  end

  setup_opts = struct('tol_rank', 1e-12);
  if isfield(cfg,'setup_opts') && ~isempty(cfg.setup_opts)
    setup_opts = cfg.setup_opts;
  end
  if ~isstruct(setup_opts)
    error('run_fetidp_block_case: cfg.setup_opts must be a struct (or empty).');
  end
  if ~isfield(setup_opts,'tol_rank'); setup_opts.tol_rank = 1e-12; end

  seed = 0;
  if isfield(cfg,'seed') && ~isempty(cfg.seed)
    seed = cfg.seed;
  end

  % Determinism (mainly to keep any randomized internal probes stable if added later)
  if exist('rng_deterministic','file') == 2
    rng_deterministic(seed);
  elseif exist('rng','file') == 2
    rng(seed);
  else
    rand('seed', seed); %#ok<RAND>
    randn('seed', seed); %#ok<RANDN>
  end

  % ------------------------------------------------------------
  % Build + setup
  % ------------------------------------------------------------
  base = build_problem_data(n, nSubX, nSubY, f_handle);
  data = setup_fetidp(base, setup_opts);

  % ------------------------------------------------------------
  % Solve dual (multiplier) problem
  % ------------------------------------------------------------
  [lambda, stats] = solve_fetidp(data, tol, maxit);

  % ------------------------------------------------------------
  % Reconstruct primal traces + global free-DOF solution
  % ------------------------------------------------------------
  [w_c, w_d, u_free, recdiag] = reconstruct_fetidp_solution(lambda, data);

  % ------------------------------------------------------------
  % Global reference (direct) and diagnostics
  % ------------------------------------------------------------
  u_ref = data.Kff \ data.Ff;

  diag = struct();
  diag.nFree      = numel(data.free);
  diag.nC         = data.primal.nC;
  diag.nDeltaProd = data.nDeltaProd;
  diag.nLambda    = data.nLambda;

  diag.rel_err_u = norm(u_free - u_ref, 2) / max(1, norm(u_ref, 2));
  diag.rel_res_global = norm(data.Kff*u_free - data.Ff, 2) / max(1, norm(data.Ff, 2));

  % Constraint residual (continuity constraints on Delta)
  if isstruct(recdiag) && isfield(recdiag,'constraint_norm')
    diag.constraint_norm = recdiag.constraint_norm;
  else
    diag.constraint_norm = norm(data.Bd * w_d, 2);
  end

  % Dual residual consistency: check A_lambda*lambda = b
  % (b is formed the same way as in solve_fetidp)
  [gC, gD] = assemble_rhs_pieces_like_solve_fetidp_(data);
  out0 = solve_tildeS(gC, gD, data);
  b = data.Bd * out0.wD;

  Al = applyA_lambda(lambda, data);
  diag.rel_res_dual = norm(Al - b, 2) / max(1, norm(b, 2));

  % Symmetry probe for A_lambda in bilinear form:
  x = ones(data.nLambda,1);
  y = (1:data.nLambda)' / max(1,data.nLambda);
  Ax = applyA_lambda(x, data);
  Ay = applyA_lambda(y, data);
  s1 = x' * Ay;
  s2 = y' * Ax;
  %diag.rel_sym_bilin = abs(s1 - s2) / max(1, abs(s1), abs(s2));
  diag.rel_sym_bilin = abs(s1 - s2) / max([1, abs(s1), abs(s2)]);

  % Positive definiteness probe (energy should be positive for nonzero x)
  diag.energy_x = x' * Ax;

  % Pack output
  out = struct();
  out.data = data;
  out.lambda = lambda;
  out.stats = stats;

  out.w_c = w_c;
  out.w_d = w_d;
  out.u_free = u_free;
  out.u_ref  = u_ref;

  out.reconstruction_diag = recdiag;
  out.diagnostics = diag;
end

% ========================================================================
% Helper: mirror RHS assembly in solve_fetidp.m
% ========================================================================
function [gC, gD] = assemble_rhs_pieces_like_solve_fetidp_(data)
  sub  = data.sub;
  nSub = numel(sub);

  gD = zeros(data.nDeltaProd, 1);
  gC = zeros(data.primal.nC, 1);

  for i = 1:nSub
    gi = sub(i).g(:);

    rng = data.delta_range{i};
    if ~isempty(rng)
      gD(rng(:)) = gi(sub(i).idx_d(:));
    end

    c_ids = sub(i).c_ids(:);
    if ~isempty(c_ids)
      gC(c_ids) = gC(c_ids) + gi(sub(i).idx_c(:));
    end
  end
end