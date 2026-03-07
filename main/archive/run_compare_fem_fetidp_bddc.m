function out = run_compare_fem_fetidp_bddc(cfg)
%RUN_COMPARE_FEM_FETIDP_BDDC  Solve same problem by FEM, FETI-DP, BDDC and compare + plot.
%
% Usage:
%   out = run_compare_fem_fetidp_bddc();                         % defaults
%   out = run_compare_fem_fetidp_bddc(struct('n',32,'nSubX',4,'nSubY',4));
%   out = run_compare_fem_fetidp_bddc(struct('f_handle',@(x,y) 2+x-3*y));
%
% Requires on path:
%   setup_paths, build_problem_data
%   setup_fetidp, solve_fetidp, reconstruct_fetidp_solution
%   setup_bddc,   solve_bddc,   reconstruct_bddc_solution

  if nargin < 1 || isempty(cfg), cfg = struct(); end
  if ~isstruct(cfg), error('cfg must be a struct or empty.'); end

  % ---------------- defaults ----------------
  if ~isfield(cfg,'n'),      cfg.n = 64; end
  if ~isfield(cfg,'nSubX'),  cfg.nSubX = 16; end
  if ~isfield(cfg,'nSubY'),  cfg.nSubY = 16; end
  if ~isfield(cfg,'seed'),   cfg.seed = 1; end
  if ~isfield(cfg,'tol'),    cfg.tol = 1e-10; end
  if ~isfield(cfg,'maxit'),  cfg.maxit = 500; end
  if ~isfield(cfg,'f_handle') || isempty(cfg.f_handle)
    cfg.f_handle = @(x,y) 1.0;
  end
  if ~isfield(cfg,'do_plots'), cfg.do_plots = true; end
  if ~isfield(cfg,'plot_diffs'), cfg.plot_diffs = false; end  % optional

  % solver-specific options (optional)
  if ~isfield(cfg,'fetidp'), cfg.fetidp = struct(); end
  if ~isfield(cfg.fetidp,'setup_opts'), cfg.fetidp.setup_opts = []; end

  if ~isfield(cfg,'bddc'), cfg.bddc = struct(); end
  if ~isfield(cfg.bddc,'setup_opts'), cfg.bddc.setup_opts = []; end

  % ---------------- paths + determinism ----------------
  if exist('setup_paths','file') ~= 2
    error('setup_paths.m not found on path.');
  end
  setup_paths();

  if exist('rng_deterministic','file') == 2
    rng_deterministic(cfg.seed);
  else
    rng(cfg.seed);
  end

  % Graphics (Octave-friendly; safe no-op in MATLAB if missing)
  try
    if exist('graphics_toolkit','file') == 2
      graphics_toolkit('qt');
    end
  catch
    % ignore
  end

  % ---------------- build shared base problem (one DD) ----------------
  base = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, cfg.f_handle);

  % FEM reference (monolithic reduced system)
  u_free_fem = base.Kff \ base.Ff;

  % ---------------- FETI-DP ----------------
  if exist('setup_fetidp','file') ~= 2 || exist('solve_fetidp','file') ~= 2 || exist('reconstruct_fetidp_solution','file') ~= 2
    error('FETI-DP functions not found (setup_fetidp / solve_fetidp / reconstruct_fetidp_solution).');
  end

  try
    if isempty(cfg.fetidp.setup_opts)
      data_f = setup_fetidp(base);
    else
      data_f = setup_fetidp(base, cfg.fetidp.setup_opts);
    end
  catch
    % fallback: some versions may only accept one arg
    data_f = setup_fetidp(base);
  end

  [lambda_f, stats_f] = solve_fetidp(data_f, cfg.tol, cfg.maxit);
  [~, ~, u_free_fetidp, diag_f] = reconstruct_fetidp_solution(lambda_f, data_f);

  % ---------------- BDDC ----------------
  if exist('setup_bddc','file') ~= 2 || exist('solve_bddc','file') ~= 2 || exist('reconstruct_bddc_solution','file') ~= 2
    error('BDDC functions not found (setup_bddc / solve_bddc / reconstruct_bddc_solution).');
  end

  try
    if isempty(cfg.bddc.setup_opts)
      data_b = setup_bddc(base);
    else
      data_b = setup_bddc(base, cfg.bddc.setup_opts);
    end
  catch
    % fallback: some versions may only accept one arg
    data_b = setup_bddc(base);
  end

  sol_b = solve_bddc(data_b, struct('tol', cfg.tol, 'maxit', cfg.maxit));
  rec_b = reconstruct_bddc_solution(sol_b.u_hat, data_b);
  u_free_bddc = rec_b.u_free;

  % ---------------- comparisons ----------------
  rel = @(a,b) norm(a-b,2) / max(1, norm(a,2));

  cmp = struct();
  cmp.rel_fetidp_vs_fem = rel(u_free_fetidp, u_free_fem);
  cmp.rel_bddc_vs_fem   = rel(u_free_bddc,   u_free_fem);
  cmp.rel_bddc_vs_fetidp= rel(u_free_bddc,   u_free_fetidp);

  cmp.maxabs_fetidp_vs_fem = max(abs(u_free_fetidp - u_free_fem));
  cmp.maxabs_bddc_vs_fem   = max(abs(u_free_bddc   - u_free_fem));
  cmp.maxabs_bddc_vs_fetidp= max(abs(u_free_bddc   - u_free_fetidp));

  fprintf('\n=== Compare solvers (n=%d, %dx%d subdomains) ===\n', cfg.n, cfg.nSubX, cfg.nSubY);
  fprintf('FEM      : direct solve on Kff\\Ff (nFree=%d)\n', numel(base.free));
  fprintf('FETI-DP  : flag=%d, iter=%d, relres=%.3e, ||Bd*wD||=%.3e\n', ...
          stats_f.flag, stats_f.iter, stats_f.relres, diag_f.constraint_norm);
  fprintf('BDDC     : flag=%d, iter=%d, relres=%.3e\n', ...
          sol_b.stats.flag, sol_b.stats.iter, sol_b.stats.relres);

  fprintf('\nErrors vs FEM:\n');
  fprintf('  rel(FETI-FEM) = %.3e,  maxabs = %.3e\n', cmp.rel_fetidp_vs_fem, cmp.maxabs_fetidp_vs_fem);
  fprintf('  rel(BDDC-FEM) = %.3e,  maxabs = %.3e\n', cmp.rel_bddc_vs_fem,   cmp.maxabs_bddc_vs_fem);
  fprintf('Consistency FETI vs BDDC:\n');
  fprintf('  rel(BDDC-FETI)= %.3e,  maxabs = %.3e\n\n', cmp.rel_bddc_vs_fetidp, cmp.maxabs_bddc_vs_fetidp);

  % ---------------- visualization ----------------
  if cfg.do_plots
    [u_fem_full, u_fet_full, u_bddc_full] = build_full_vectors_(base, u_free_fem, u_free_fetidp, u_free_bddc);

    % shared color scale
    zmin = min([u_fem_full; u_fet_full; u_bddc_full]);
    zmax = max([u_fem_full; u_fet_full; u_bddc_full]);

    figure('Name','Solutions: FEM vs FETI-DP vs BDDC');
    subplot(1,3,1); plot_solution_(base, u_fem_full,  sprintf('FEM (monolithic)'));
    caxis([zmin zmax]);

    subplot(1,3,2); plot_solution_(base, u_fet_full,  sprintf('FETI-DP (iter=%d)', stats_f.iter));
    caxis([zmin zmax]);

    subplot(1,3,3); plot_solution_(base, u_bddc_full, sprintf('BDDC (iter=%d)', sol_b.stats.iter));
    caxis([zmin zmax]);

    if cfg.plot_diffs
      figure('Name','Differences (full nodal vectors)');
      subplot(1,2,1); plot_solution_(base, u_fet_full - u_fem_full,  'FETI-DP - FEM'); colorbar;
      subplot(1,2,2); plot_solution_(base, u_bddc_full - u_fem_full, 'BDDC - FEM');   colorbar;
    end
  end

  % ---------------- pack output ----------------
  out = struct();
  out.cfg = cfg;
  out.base = base;

  out.fem.u_free = u_free_fem;

  out.fetidp.lambda = lambda_f;
  out.fetidp.stats  = stats_f;
  out.fetidp.diag   = diag_f;
  out.fetidp.u_free = u_free_fetidp;

  out.bddc.sol   = sol_b;
  out.bddc.rec   = rec_b;
  out.bddc.u_free= u_free_bddc;

  out.compare = cmp;
end

% =======================================================================
% local helpers
% =======================================================================
function [u_fem, u_fet, u_bddc] = build_full_vectors_(base, u_free_fem, u_free_fet, u_free_bddc)
  N = size(base.p,1);
  u_fem  = zeros(N,1); u_fem(base.free)  = u_free_fem;
  u_fet  = zeros(N,1); u_fet(base.free)  = u_free_fet;
  u_bddc = zeros(N,1); u_bddc(base.free) = u_free_bddc;
end

function plot_solution_(base, u, ttl)
  trisurf(base.t, base.p(:,1), base.p(:,2), u);
  view(3); grid on; axis vis3d;
  shading interp; colorbar;
  xlabel('x'); ylabel('y'); zlabel('u_h');
  title(ttl);
end