function out = run_fem_block_case(cfg)
%RUN_FEM_BLOCK_CASE  End-to-end FEM block run (mesh -> assemble -> eliminate -> solve) with diagnostics.
%
% This helper is intended for intermediate/integration testing of the FEM block.
%
% Required fields:
%   cfg.n        - mesh resolution (integer >= 1)
%   cfg.f_handle - RHS handle f(x,y) for -Delta u = f
%
% Optional fields:
%   cfg.mutate          - function handle [p,t,bnd] = mutate(p,t,bnd)
%                         called after mesh generation (fault injection).
%   cfg.u_exact_handle  - function handle u(x,y) returning exact solution at points.
%   cfg.solve           - logical (default true). If false, stop after elimination.
%
% Output struct out contains:
%   p,t,bnd,K,F,Kff,Ff,free_nodes,u_full,u_free and diagnostics.

  if nargin ~= 1
    error('run_fem_block_case: expected one input struct cfg.');
  end
  if ~isstruct(cfg)
    error('run_fem_block_case: cfg must be a struct.');
  end

  if ~isfield(cfg,'n')
    error('run_fem_block_case: cfg.n is required.');
  end
  if ~isfield(cfg,'f_handle')
    error('run_fem_block_case: cfg.f_handle is required.');
  end

  n = cfg.n;
  f_handle = cfg.f_handle;

  do_solve = true;
  if isfield(cfg,'solve')
    do_solve = logical(cfg.solve);
  end

  % Mesh
  [p, t, bnd] = mesh_unit_square_P1(n);

  % Optional mutation hook (fault injection / degenerate cases)
  if isfield(cfg,'mutate') && ~isempty(cfg.mutate)
    if ~isa(cfg.mutate,'function_handle')
      error('run_fem_block_case: cfg.mutate must be a function handle.');
    end
    [p, t, bnd] = cfg.mutate(p, t, bnd);
  end

  % Assembly
  K = assemble_stiffness_P1(p, t);
  F = assemble_load_P1(p, t, f_handle);

  % Dirichlet elimination
  if ~isfield(bnd,'dirichlet_nodes')
    error('run_fem_block_case: bnd.dirichlet_nodes missing.');
  end
  [Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, F, bnd.dirichlet_nodes);

  u_free = [];
  u_full = [];
  if do_solve
    u_free = Kff \ Ff;
    u_full = zeros(size(p,1), 1);
    u_full(free_nodes) = u_free;
  end

  % Diagnostics
  diag = struct();

  diag.N  = size(p,1);
  diag.Ne = size(t,1);
  diag.nFree = numel(free_nodes);

  % Symmetry residuals (scale-aware)
  diag.rel_sym_K   = norm(K - K.', 'fro') / max(1, norm(K, 'fro'));
  diag.rel_sym_Kff = norm(Kff - Kff.', 'fro') / max(1, norm(Kff, 'fro'));

  % SPD check (Cholesky)
  diag.chol_ok = false;
  try
    chol(Kff);
    diag.chol_ok = true;
  catch
    diag.chol_ok = false;
  end

  % Constant nullspace check on unreduced K
  one = ones(size(K,1),1);
  diag.rel_null_K = norm(K*one,2) / max(1, norm(K,1) * norm(one,2));

  % Solve residual
  if do_solve
    r = Kff*u_free - Ff;
    diag.rel_resid = norm(r,2) / max(1, norm(Ff,2));
  else
    diag.rel_resid = NaN;
  end

  % Simple probe values (bottom/top center nodes if they exist)
  diag.probe = struct();
  np = bnd.n + 1;
  if mod(bnd.n,2) == 0
    i_mid = bnd.n/2 + 1;
    idx_bot = i_mid;                  % j=1
    idx_top = i_mid + bnd.n*np;       % j=np
    diag.probe.idx_bot = idx_bot;
    diag.probe.idx_top = idx_top;
    diag.probe.xy_bot  = p(idx_bot,:);
    diag.probe.xy_top  = p(idx_top,:);
    if do_solve
      diag.probe.u_bot = u_full(idx_bot);
      diag.probe.u_top = u_full(idx_top);
    else
      diag.probe.u_bot = NaN;
      diag.probe.u_top = NaN;
    end
  else
    diag.probe.idx_bot = NaN;
    diag.probe.idx_top = NaN;
  end

  % Exact-solution errors (if provided)
  if isfield(cfg,'u_exact_handle') && ~isempty(cfg.u_exact_handle)
    if ~isa(cfg.u_exact_handle,'function_handle')
      error('run_fem_block_case: cfg.u_exact_handle must be a function handle.');
    end
    if do_solve
      uex = cfg.u_exact_handle(p(:,1), p(:,2));
      if ~isvector(uex) || numel(uex) ~= size(p,1)
        error('run_fem_block_case: u_exact_handle must return N values for N nodes.');
      end
      uex = uex(:);

      % Errors on free nodes (most informative)
      ef = u_full(free_nodes) - uex(free_nodes);
      uf = uex(free_nodes);

      diag.errL2_abs = norm(ef,2);
      diag.errL2_rel = diag.errL2_abs / max(1e-16, norm(uf,2));

      diag.errInf_abs = norm(ef,inf);
      diag.errInf_rel = diag.errInf_abs / max(1e-16, norm(uf,inf));
    else
      diag.errL2_abs = NaN;
      diag.errL2_rel = NaN;
      diag.errInf_abs = NaN;
      diag.errInf_rel = NaN;
    end
  end

  % Pack output
  out = struct();
  out.p = p;
  out.t = t;
  out.bnd = bnd;
  out.K = K;
  out.F = F;
  out.Kff = Kff;
  out.Ff = Ff;
  out.free_nodes = free_nodes;
  out.u_free = u_free;
  out.u_full = u_full;
  out.diagnostics = diag;
end