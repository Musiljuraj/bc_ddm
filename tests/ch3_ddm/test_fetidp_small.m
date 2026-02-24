% ============================================================
% File: tests/ch4_feti_dp/test_fetidp_small.m
% ============================================================
function test_fetidp_small()
%TEST_FETIDP_SMALL  Deterministic smoke test for FETI-DP end-to-end.
%
% Checks:
%   - PCG converges (flag==0)
%   - constraint residual small
%   - reconstructed u_free close to monolithic u_ref (Kff\Ff)

  here = fileparts(mfilename('fullpath'));          % .../bc_ddm/tests/ch4_feti_dp
  root = fullfile(here, '..', '..');                % .../bc_ddm
  addpath(fullfile(root, 'main'));                  % make setup_paths visible
  setup_paths();                                    % now adds src/ and tests/
  rng_deterministic(1);

  n = 12;
  nSubX = 3;
  nSubY = 3;
  f_handle = @(x,y) 1.0;

  data = build_problem_data(n, nSubX, nSubY, f_handle);
  data = setup_fetidp(data, struct('tol_rank', 1e-12));

  [lambda, stats] = solve_fetidp(data, 1e-8, 500);
  [~, ~, u_free, diag] = reconstruct_fetidp_solution(lambda, data);

  u_ref = data.Kff \ data.Ff;

  relerr = norm(u_free - u_ref) / max(1.0, norm(u_ref));

  assert(stats.flag == 0, 'FETI-DP PCG did not converge (flag=%d).', stats.flag);
  assert(diag.constraint_norm < 1e-6, 'Constraint residual too large: %.3e', diag.constraint_norm);
  assert(relerr < 1e-4, 'Relative error too large: %.3e', relerr);

  fprintf('PASS: test_fetidp_small (iter=%d, relerr=%.3e, ||Bd*wD||=%.3e)\n', ...
          stats.iter, relerr, diag.constraint_norm);
end