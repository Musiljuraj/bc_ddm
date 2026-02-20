% ============================================================
% File: main/run_fetidp_demo.m
% ============================================================
function results = run_fetidp_demo()
%RUN_FETIDP_DEMO  Minimal end-to-end FETI-DP demo on the unit square.
%
% This demo:
%   - builds data via build_problem_data
%   - sets up FETI-DP via setup_fetidp
%   - solves for lambda with PCG in multiplier space
%   - reconstructs u_free and compares to monolithic reference u_ref = Kff\Ff
%
% Adjust (n, nSubX, nSubY) here.

  setup_paths();

  % Parameters (default consistent with structured partitioning):
  n = 16;        % mesh intervals per axis => (n+1)x(n+1) nodes
  nSubX = 4;
  nSubY = 4;

  f_handle = @(x,y) -1.0;

  data = build_problem_data(n, nSubX, nSubY, f_handle);
  data = setup_fetidp(data, struct('tol_rank', 1e-12));

  tol = 1e-8;
  maxit = 500;

  [lambda, stats] = solve_fetidp(data, tol, maxit);
  [w_c, w_d, u_free, diag] = reconstruct_fetidp_solution(lambda, data); %#ok<ASGLU>

  u_ref = data.Kff \ data.Ff;

  relerr = norm(u_free - u_ref) / max(1.0, norm(u_ref));

  results = struct();
  results.stats = stats;
  results.relerr = relerr;
  results.constraint_norm = diag.constraint_norm;

  fprintf('FETI-DP demo:\n');
  fprintf('  flag=%d, iter=%d, relres=%.3e\n', stats.flag, stats.iter, stats.relres);
  fprintf('  relerr(u_free vs u_ref) = %.3e\n', relerr);
  fprintf('  ||Bd*wD||_2 = %.3e\n', diag.constraint_norm);
end