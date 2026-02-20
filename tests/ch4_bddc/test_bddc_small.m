% ============================================================
% File: tests/ch4_bddc/test_bddc_small.m
% ============================================================
function test_bddc_small()
%TEST_BDDC_SMALL  Deterministic end-to-end test of baseline BDDC.

  % Ensure project paths.
  here = fileparts(mfilename('fullpath'));
  root = fullfile(here, '..', '..');
  addpath(fullfile(root, 'main'));
  setup_paths();

  rng_deterministic(1);

  % Small deterministic instance.
  n = 8;
  nSubX = 2;
  nSubY = 2;

  data = build_problem_data(n, nSubX, nSubY, []);
  data = setup_bddc(data);

  % Symmetry sanity for A_hat.
  x = rand(data.prod.nHat, 1);
  y = rand(data.prod.nHat, 1);
  Ax = applyA_hat(x, data);
  Ay = applyA_hat(y, data);
  sym_err = abs(x' * Ay - y' * Ax) / max(1, abs(x' * Ay));
  assert(sym_err < 1e-10, 'A_hat symmetry check failed.');

  % Solve.
  opts = struct('tol', 1e-10, 'maxit', 400);
  sol = solve_bddc(data, opts);
  assert(sol.stats.flag == 0, 'BDDC PCG did not converge (flag=%d).', sol.stats.flag);

  % Reconstruct and compare to monolithic reference.
  rec = reconstruct_bddc_solution(sol.u_hat, data);
  u_ref = data.Kff \ data.Ff;
  relerr = norm(rec.u_free - u_ref) / norm(u_ref);
  assert(relerr < 1e-8, 'BDDC solution error too large: %.3e', relerr);

  fprintf('PASS: test_bddc_small (iter=%d, relerr=%.3e)\n', sol.stats.iter, relerr);
end

% Allow running as a script.
if ~isdeployed && isempty(getenv('OCTAVE_TEST_DISABLE_AUTO'))
  try
    test_bddc_small();
  catch err
    rethrow(err);
  end
end