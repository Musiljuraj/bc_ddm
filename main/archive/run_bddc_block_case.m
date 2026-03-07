function out = run_bddc_block_case(cfg)
%RUN_BDDC_BLOCK_CASE  End-to-end BDDC block run with strong diagnostics.
%
% Chain:
%   build_problem_data -> setup_bddc -> (explicit A_hat build) ->
%   applyA_hat/applyM_bddc consistency -> solve_bddc -> reconstruct_bddc_solution
%
% Required fields:
%   cfg.n, cfg.nSubX, cfg.nSubY
%
% Optional:
%   cfg.f_handle     (default [])
%   cfg.seed         (default 1)
%   cfg.opts_bddc    (default struct(); e.g. struct('tol_chol',0))
%   cfg.opts_solve   (default struct('tol',1e-10,'maxit',400))
%   cfg.mutate_data  function handle data = mutate_data(data) (fault injection)
%   cfg.solve        logical (default true)
%
% Output:
%   out.data, out.bddc_explicit, out.sol, out.rec, out.diagnostics

  if nargin ~= 1 || ~isstruct(cfg)
    error('run_bddc_block_case: expected one input struct cfg.');
  end

  req = {'n','nSubX','nSubY'};
  for k = 1:numel(req)
    if ~isfield(cfg, req{k})
      error('run_bddc_block_case: cfg.%s is required.', req{k});
    end
  end

  % Defaults
  if ~isfield(cfg,'seed') || isempty(cfg.seed), cfg.seed = 1; end
  if exist('rng_deterministic','file') == 2
    rng_deterministic(cfg.seed);
  end

  if ~isfield(cfg,'f_handle'), cfg.f_handle = []; end

  do_solve = true;
  if isfield(cfg,'solve') && ~isempty(cfg.solve)
    do_solve = logical(cfg.solve);
  end

  opts_bddc = struct();
  if isfield(cfg,'opts_bddc') && ~isempty(cfg.opts_bddc)
    opts_bddc = cfg.opts_bddc;
  end

  opts_solve = struct('tol', 1e-10, 'maxit', 400);
  if isfield(cfg,'opts_solve') && ~isempty(cfg.opts_solve)
    opts_solve = cfg.opts_solve;
    if ~isfield(opts_solve,'tol'),   opts_solve.tol = 1e-10; end
    if ~isfield(opts_solve,'maxit'), opts_solve.maxit = 400; end
  end

  % ------------------------------------------------------------
  % Build data (Ch2+Ch3 pipeline wrapper)
  % ------------------------------------------------------------
  data = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, cfg.f_handle);

  % Fault injection hook (mutate data BEFORE setup_bddc)
  if isfield(cfg,'mutate_data') && ~isempty(cfg.mutate_data)
    if ~isa(cfg.mutate_data,'function_handle')
      error('run_bddc_block_case: cfg.mutate_data must be a function handle.');
    end
    data = cfg.mutate_data(data);
  end

  % ------------------------------------------------------------
  % Setup BDDC (Ch4)
  % ------------------------------------------------------------
  data = setup_bddc(data, opts_bddc);

  % ------------------------------------------------------------
  % Explicit operator build in product/hat spaces (for strong checks)
  % ------------------------------------------------------------
  R = data.bddc.R;
  g_prod = data.bddc.g_prod(:);

  nProd = size(R,1);
  nHat  = size(R,2);

  % Assemble explicit S_prod = diag(S_i) using assembled local Schurs sub(i).S
  S_prod = sparse(nProd, nProd);
  sub = data.sub;
  for i = 1:numel(sub)
    idx = sub(i).prod_idx(:);
    if isempty(idx), continue; end
    if ~isfield(sub(i),'S') || isempty(sub(i).S)
      error('run_bddc_block_case: sub(%d).S missing/empty (need assemble_S=true upstream).', i);
    end
    S_prod(idx, idx) = sub(i).S;
  end

  A_hat = R' * (S_prod * R);
  f_hat = R' * g_prod;

  % ------------------------------------------------------------
  % Diagnostics (matrix-free vs explicit identities + symmetry checks)
  % ------------------------------------------------------------
  d = struct();
  d.nHat  = nHat;
  d.nProd = nProd;
  d.nC    = data.primal.nC;
  d.Nsub  = numel(data.sub);

  if nHat == 0
    d.rel_sym_A   = 0;
    d.rel_applyA  = 0;
    d.rel_sym_M   = 0;
    d.M_pos_quad  = 0;
  else
    d.rel_sym_A = norm(A_hat - A_hat.', 'fro') / max(1, norm(A_hat, 'fro'));

    x = (1:nHat).' / max(1, nHat);
    y_mf  = applyA_hat(x, data);
    y_ref = A_hat * x;
    d.rel_applyA = norm(y_mf - y_ref, 2) / max(1, norm(y_ref, 2));

    % Preconditioner symmetry sanity in bilinear form:
    %   x' M y  ==  y' M x  (approx)
    if exist('rng_deterministic','file') == 2
      rng_deterministic(12345);
    end
    x = rand(nHat, 1);
    y = rand(nHat, 1);
    Mx = applyM_bddc(x, data);
    My = applyM_bddc(y, data);
    d.rel_sym_M = abs(x' * My - y' * Mx) / max(1, abs(x' * My));

    d.M_pos_quad = x' * Mx; % should be > 0 for nonzero x if SPD
  end

  % ------------------------------------------------------------
  % Solve + reconstruct + compare to monolithic reference
  % ------------------------------------------------------------
  sol = [];
  rec = [];
  if do_solve
    sol = solve_bddc(data, opts_solve);

    % PCG residual w.r.t. explicit A_hat (strong, deterministic)
    if nHat == 0
      d.rel_res_pcg = 0;
      d.rel_u_hat_direct = 0;
    else
      r = A_hat * sol.u_hat - f_hat;
      d.rel_res_pcg = norm(r, 2) / max(1, norm(f_hat, 2));

      % Direct reference in hat-space (only for moderate sizes)
      if nHat <= 400
        u_dir = A_hat \ f_hat;
        d.rel_u_hat_direct = norm(sol.u_hat - u_dir, 2) / max(1, norm(u_dir, 2));
      else
        d.rel_u_hat_direct = NaN;
      end
    end

    rec = reconstruct_bddc_solution(sol.u_hat, data);
    u_ref = data.Kff \ data.Ff;
    d.rel_u_free_ref = norm(rec.u_free - u_ref, 2) / max(1, norm(u_ref, 2));
  end

  out = struct();
  out.data = data;

  out.bddc_explicit = struct();
  out.bddc_explicit.S_prod = S_prod;
  out.bddc_explicit.A_hat  = A_hat;
  out.bddc_explicit.f_hat  = f_hat;

  out.sol = sol;
  out.rec = rec;
  out.diagnostics = d;
end