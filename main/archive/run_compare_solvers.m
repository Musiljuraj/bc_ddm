% ============================================================
% File: main/run_compare_solvers.m
% ============================================================
function out = run_compare_solvers()
%RUN_COMPARE_SOLVERS  Compare FETI-DP and BDDC on the same deterministic instance.
%
% If BDDC is not available, this will still run FETI-DP and the monolithic
% reference solve.

  setup_paths();

  % Shared deterministic configuration.
  n = 16;
  nSubX = 4;
  nSubY = 4;
  seed = 1;
  rng_deterministic(seed);

  data = build_problem_data(n, nSubX, nSubY, []);
  u_ref = data.Kff \ data.Ff;

  % --- FETI-DP ---
  data_f = setup_fetidp(data);

  
  tol_f = 1e-8;
  maxit_f = 500;
  [lambda_f, stats_f] = solve_fetidp(data_f, tol_f, maxit_f);

  sol_f = struct();
  sol_f.lambda = lambda_f;
  sol_f.stats  = stats_f;
  %opts_f = struct('tol', 1e-8, 'maxit', 500);
  %sol_f = solve_fetidp(data_f, opts_f.tol, opts_f.maxit);
  
  [~, ~, u_free_f, diag_f] = reconstruct_fetidp_solution(sol_f.lambda, data_f);
  relerr_f = norm(u_free_f - u_ref) / norm(u_ref);  
  rec_f = struct();
  rec_f.u_free = u_free_f;
  rec_f.diag   = diag_f;
  %rec_f = reconstruct_fetidp_solution(sol_f.lambda, data_f);
  %relerr_f = norm(rec_f.u_free - u_ref) / norm(u_ref);


  fprintf('FETI-DP: flag=%d, iter=%d, relres=%.3e, relerr(u)=%.3e, ||Bd*wD||=%.3e\n', ...
        sol_f.stats.flag, sol_f.stats.iter, sol_f.stats.relres, relerr_f, diag_f.constraint_norm);
  %fprintf('FETI-DP: flag=%d, iter=%d, relres=%.3e, relerr(u)=%.3e\n', ...
          %sol_f.stats.flag, sol_f.stats.iter, sol_f.stats.relres, relerr_f);

  % --- BDDC (if present) ---
  have_bddc = exist('setup_bddc','file') == 2 && exist('solve_bddc','file') == 2;

  if have_bddc
    data_b = setup_bddc(data);
    opts_b = struct('tol', 1e-8, 'maxit', 500);
    sol_b = solve_bddc(data_b, opts_b);
    rec_b = reconstruct_bddc_solution(sol_b.u_hat, data_b);
    relerr_b = norm(rec_b.u_free - u_ref) / norm(u_ref);

    fprintf('BDDC   : flag=%d, iter=%d, relres=%.3e, relerr(u)=%.3e\n', ...
            sol_b.stats.flag, sol_b.stats.iter, sol_b.stats.relres, relerr_b);
  else
    sol_b = [];
    rec_b = [];
    relerr_b = NaN;
    fprintf('BDDC: not available in path (skipped).\n');
  end

  out = struct();
  out.n = n;
  out.nSubX = nSubX;
  out.nSubY = nSubY;
  out.seed = seed;
  out.u_ref = u_ref;
  out.fetidp.sol = sol_f;
  out.fetidp.rec = rec_f;
  out.fetidp.relerr = relerr_f;
  out.bddc.sol = sol_b;
  out.bddc.rec = rec_b;
  out.bddc.relerr = relerr_b;
end