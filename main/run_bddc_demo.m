% ============================================================
% File: main/run_bddc_demo.m
% ============================================================
function out = run_bddc_demo()
%RUN_BDDC_DEMO  Minimal end-to-end demo of baseline BDDC on the model problem.
%
% This demo builds the problem data, sets up BDDC, solves the hat-space system
% by PCG, reconstructs the global solution, and compares to a monolithic
% reference solve on Kff.

  setup_paths();

  % Deterministic configuration (edit as needed).
  cfg = struct();
  cfg.n = 16;        % mesh parameter (elements per side)
  cfg.nSubX = 4;
  cfg.nSubY = 4;
  cfg.seed = 1;

  rng_deterministic(cfg.seed);

  % Build reusable problem data (Ch.2 + Ch.3 pipeline).
  data = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, []);

  % Setup BDDC.
  data = setup_bddc(data);

  % Solve.
  opts = struct('tol', 1e-8, 'maxit', 500);
  sol = solve_bddc(data, opts);

  % Reconstruct global solution.
  rec = reconstruct_bddc_solution(sol.u_hat, data);

  % Reference monolithic solve.
  u_ref = data.Kff \ data.Ff;
  relerr = norm(rec.u_free - u_ref) / norm(u_ref);

  fprintf('BDDC PCG: flag=%d, iter=%d, relres=%.3e, relerr(u)=%.3e\n', ...
          sol.stats.flag, sol.stats.iter, sol.stats.relres, relerr);

  out = struct();
  out.cfg = cfg;
  out.sol = sol;
  out.rec = rec;
  out.u_ref = u_ref;
  out.relerr = relerr;
end