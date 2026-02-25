function test_solve_bddc()
%TEST_SOLVE_BDDC  Basic unit test for solve_bddc (smoke + consistency).
%
% Tests the function: solve_bddc
% This test verifies:
%   - wrong-arg-count throws
%   - missing data.bddc throws with expected message
%   - can build a tiny problem and run solve_bddc deterministically
%   - output sizes are consistent with R; outputs are finite
%   - consistency: applyA_hat(u_hat) approximately equals f_hat = R'*g_prod
%   - consistency: w_prod equals R*u_hat

  ensure_project_paths_();

  fprintf('RUNNING: test_solve_bddc (tests solve_bddc)\n');
  fprintf('  Checks: arg-count error; missing bddc error; output size/finite; residual applyA_hat(u)=f_hat; w_prod=R*u.\n');

  % -----------------------
  % Negative: wrong arg count
  % -----------------------
  assert_throws_(@() solve_bddc());

  % -----------------------
  % Negative: missing data.bddc
  % -----------------------
  assert_throws_msg_( ...
    @() solve_bddc(struct(), struct()), ...
    'data.bddc missing' );

  % -----------------------
  % Build tiny BDDC data
  % -----------------------
  data = make_tiny_bddc_data_();

  % Basic structural sanity (minimal, cheap)
  assert_true_(isstruct(data), 'data must be a struct.');
  assert_true_(isfield(data,'bddc') && isstruct(data.bddc), 'Missing/invalid data.bddc.');
  assert_true_(isfield(data.bddc,'R') && isnumeric(data.bddc.R), 'Missing/invalid data.bddc.R.');

  R = data.bddc.R;
  assert_true_(ismatrix(R) && size(R,1) > 0 && size(R,2) > 0, 'data.bddc.R must be nonempty 2D.');
  nProd = size(R,1);
  nHat  = size(R,2);

  % -----------------------
  % Deterministic RHS in product space (override g_prod only)
  % -----------------------
  data.bddc.g_prod = (1:nProd)';  % deterministic, nontrivial
  g_prod = data.bddc.g_prod;

  % Form f_hat exactly as solve_bddc does
  f_hat = R' * g_prod;
  assert_true_(isvector(f_hat) && numel(f_hat) == nHat, 'f_hat size mismatch vs size(R,2).');
  assert_true_(isreal(f_hat) && all(isfinite(f_hat)), 'f_hat must be real and finite.');

  % -----------------------
  % Run solver
  % -----------------------
  opts = struct('tol', 1e-10, 'maxit', 300);
  out = solve_bddc(data, opts);

  % Output checks
  assert_true_(isstruct(out), 'solve_bddc must return a struct.');
  assert_true_(isfield(out,'u_hat') && isnumeric(out.u_hat), 'out.u_hat missing/invalid.');
  assert_true_(isfield(out,'w_prod') && isnumeric(out.w_prod), 'out.w_prod missing/invalid.');
  assert_true_(isfield(out,'stats'), 'out.stats missing.');

  u_hat = out.u_hat(:);
  w_prod = out.w_prod(:);

  assert_true_(isvector(u_hat) && numel(u_hat) == nHat, 'u_hat length mismatch vs size(R,2).');
  assert_true_(isreal(u_hat) && all(isfinite(u_hat)), 'u_hat must be real and finite.');
  assert_true_(isvector(w_prod) && numel(w_prod) == nProd, 'w_prod length mismatch vs size(R,1).');
  assert_true_(isreal(w_prod) && all(isfinite(w_prod)), 'w_prod must be real and finite.');

  % stats.flag == 0 when available (do not require specific fields)
  if isstruct(out.stats) && isfield(out.stats,'flag')
    assert_true_(out.stats.flag == 0, 'PCG did not report convergence (stats.flag ~= 0).');
  end

  % -----------------------
  % Consistency: w_prod = R*u_hat
  % -----------------------
  w_ref = R * u_hat;
  assert_true_(numel(w_ref) == numel(w_prod), 'R*u_hat size mismatch vs w_prod.');
  rel_w = norm(w_prod - w_ref) / max(1, norm(w_ref));
  assert_true_(rel_w <= 1e-12, 'w_prod not consistent with R*u_hat.');

  % -----------------------
  % Consistency: applyA_hat(u_hat) ~= f_hat
  % -----------------------
  assert_true_(exist('applyA_hat','file') == 2, 'applyA_hat not found on path.');
  Au = applyA_hat(u_hat, data);
  assert_true_(isvector(Au) && numel(Au) == nHat, 'applyA_hat output size mismatch.');
  assert_true_(isreal(Au) && all(isfinite(Au)), 'applyA_hat output must be real and finite.');

  rel_res = norm(Au - f_hat) / max(1, norm(f_hat));
  assert_true_(rel_res <= 5e-7, 'Relative residual too large for applyA_hat(u_hat) = f_hat.');

  fprintf('PASS: test_solve_bddc — tested solve_bddc\n');
  fprintf('  Verified: arg-count throws; missing bddc throws; u_hat/w_prod size+finite; relres=%.3e; w_rel=%.3e\n', rel_res, rel_w);

end

% =====================================================================
% Local helpers
% =====================================================================

function data = make_tiny_bddc_data_()
  assert_true_(exist('build_problem_data','file') == 2, 'build_problem_data not found on path.');
  assert_true_(exist('setup_bddc','file') == 2, 'setup_bddc not found on path.');

  % Try a few small structured decompositions; stop at first that works.
  cases = {
    [2, 2, 1]
    [2, 1, 2]
    [2, 2, 2]
    [4, 2, 2]
    [6, 2, 1]
    [6, 3, 1]
  };

  lastErr = [];
  for k = 1:numel(cases)
    c = cases{k};
    n     = c(1);
    nSubX = c(2);
    nSubY = c(3);

    try
      base = build_problem_data(n, nSubX, nSubY);

      try
        data = setup_bddc(base);
      catch
        data = setup_bddc(base, struct());
      end

      % Minimal requirements for solve_bddc pipeline
      if isstruct(data) && isfield(data,'bddc') && isstruct(data.bddc) ...
            && isfield(data.bddc,'R') && ~isempty(data.bddc.R)
        return;
      end
    catch err
      lastErr = err;
    end
  end

  if ~isempty(lastErr), rethrow(lastErr); end
  error('test_solve_bddc:setup_failed', 'Could not build a suitable tiny BDDC configuration.');
end

function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did
    return;
  end

  % Strategy 1: setup_paths already resolvable
  if exist('setup_paths','file') == 2
    setup_paths();
    sp = which('setup_paths');
    main_dir = fileparts(sp);
    root_dir = fileparts(main_dir);
    addpath(genpath(fullfile(root_dir,'tests')));
    did = true;
    return;
  end

  % Strategy 2: derive root from this test file location
  thisdir = fileparts(mfilename('fullpath'));
  root_dir = thisdir;

  for k = 1:10
    cand = fullfile(root_dir, 'main', 'setup_paths.m');
    if exist(cand, 'file') == 2
      addpath(fullfile(root_dir,'main'));
      setup_paths();
      addpath(genpath(fullfile(root_dir,'tests')));
      did = true;
      return;
    end
    parent = fileparts(root_dir);
    if isempty(parent) || strcmp(parent, root_dir)
      break;
    end
    root_dir = parent;
  end

  error('test_solve_bddc:paths', 'could not locate project root containing main/setup_paths.m');
end

function assert_true_(cond, msg, varargin)
  if nargin < 2, msg = 'assertion failed'; end
  if ~cond
    if nargin > 2
      msg = sprintf(msg, varargin{:});
    end
    error('test_solve_bddc:assert', '%s', msg);
  end
end

function assert_throws_(fh)
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  if ~threw
    error('test_solve_bddc:expected_error', 'Expected an error but none was thrown.');
  end
end

function assert_throws_msg_(fh, msg_substr)
  threw = false;
  got_msg = '';
  try
    fh();
  catch err
    threw = true;
    got_msg = err.message;
  end
  if ~threw
    error('test_solve_bddc:expected_error', 'Expected an error but none was thrown.');
  end
  if nargin >= 2 && ~isempty(msg_substr)
    assert_true_(~isempty(strfind(got_msg, msg_substr)), ...
      'Error message did not contain expected substring "%s". Got: %s', msg_substr, got_msg);
  end
end