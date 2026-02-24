function test_solve_fetidp()
%TEST_SOLVE_FETIDP  Basic unit test for solve_fetidp (smoke + consistency).
%
% Tests the function: solve_fetidp
% This test verifies:
%   - wrong-arg-count throws
%   - can build a tiny problem and run solve_fetidp deterministically
%   - output has correct length and finite values
%   - lambda approximately satisfies A_lambda*lambda = b, where b is formed
%     in the same way as solve_fetidp (via solve_tildeS + Bd*wD0)

  ensure_project_paths_();

  fprintf('RUNNING: test_solve_fetidp (tests solve_fetidp)\n');
  fprintf('  Checks: arg-count error; output size/finite; residual A_lambda*lambda=b.\n');

  % -----------------------
  % Negative: wrong arg count
  % -----------------------
  assert_throws_(@() solve_fetidp());

  % -----------------------
  % Build tiny FETI-DP data
  % -----------------------
  data = make_tiny_fetidp_data_();

  % Basic structural sanity
  assert_true_(isstruct(data), 'data must be a struct.');
  assert_true_(isfield(data,'sub') && numel(data.sub) >= 2, 'Need >=2 subdomains.');
  assert_true_(isfield(data,'primal') && isfield(data.primal,'nC'), 'Missing data.primal.nC.');
  assert_true_(isfield(data,'nLambda') && data.nLambda > 0, 'Missing/invalid data.nLambda.');
  assert_true_(isfield(data,'nDeltaProd') && data.nDeltaProd >= 0, 'Missing/invalid data.nDeltaProd.');
  assert_true_(isfield(data,'Bd') && isnumeric(data.Bd), 'Missing/invalid data.Bd.');
  assert_true_(isfield(data,'BdT') && isnumeric(data.BdT), 'Missing/invalid data.BdT.');

  % -----------------------
  % Force deterministic RHS sub(i).g
  % -----------------------
  data = set_deterministic_local_rhs_(data);

  % Assemble RHS pieces exactly like solve_fetidp
  [gC, gD] = assemble_rhs_pieces_like_solve_fetidp_(data);

  % Form expected multiplier RHS b the same way solve_fetidp does
  out0 = solve_tildeS(gC, gD, data);
  assert_true_(isstruct(out0) && isfield(out0,'wD'), 'solve_tildeS must return struct with field wD.');
  b = data.Bd * out0.wD;

  assert_true_(isvector(b) && numel(b) == data.nLambda, 'b must be a vector of length data.nLambda.');
  assert_true_(isreal(b) && all(isfinite(b)), 'b must be real and finite.');

  % -----------------------
  % Run solver
  % -----------------------
  tol = 1e-10;
  maxit = 300;
  [lambda, stats] = solve_fetidp(data, tol, maxit);

  % Output checks
  assert_true_(isnumeric(lambda) && isvector(lambda), 'lambda must be a numeric vector.');
  lambda = lambda(:);
  assert_true_(numel(lambda) == data.nLambda, 'lambda length mismatch vs data.nLambda.');
  assert_true_(isreal(lambda) && all(isfinite(lambda)), 'lambda must be real and finite.');

  if isstruct(stats)
    if isfield(stats,'flag')
      assert_true_(stats.flag == 0, 'PCG did not report convergence (stats.flag ~= 0).');
    end
  end

  % -----------------------
  % Consistency: residual in multiplier space
  % -----------------------
  Al = applyA_lambda(lambda, data);
  assert_true_(isvector(Al) && numel(Al) == data.nLambda, 'applyA_lambda output size mismatch.');
  assert_true_(isreal(Al) && all(isfinite(Al)), 'applyA_lambda output must be real and finite.');

  rel_res = norm(Al - b) / max(1, norm(b));
  assert_true_(rel_res <= 5e-7, 'Relative residual too large for A_lambda*lambda = b.');

  fprintf('PASS: test_solve_fetidp — tested solve_fetidp\n');
  fprintf('  Verified: arg-count throws; lambda size/finite; relres(A_lambda*lambda,b)=%.3e\n', rel_res);

end

% =====================================================================
% Local helpers
% =====================================================================

function data = make_tiny_fetidp_data_()
  assert_true_(exist('build_problem_data','file') == 2, 'build_problem_data not found on path.');
  assert_true_(exist('setup_fetidp','file') == 2, 'setup_fetidp not found on path.');

  % build_problem_data(n, nSubX, nSubY, f_handle) with narginchk(3,4)
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
        data = setup_fetidp(base);
      catch
        data = setup_fetidp(base, struct());
      end

      % Minimal requirements for solve_fetidp pipeline
      if isstruct(data) && isfield(data,'sub') && numel(data.sub) >= 2 ...
            && isfield(data,'primal') && isfield(data.primal,'nC') ...
            && isfield(data,'Bd') && isfield(data,'BdT') ...
            && isfield(data,'nLambda') && isfield(data,'nDeltaProd') ...
            && data.nLambda > 0
        return;
      end
    catch err
      lastErr = err;
    end
  end

  if ~isempty(lastErr), rethrow(lastErr); end
  error('test_solve_fetidp:setup_failed', 'Could not build a suitable tiny FETI-DP configuration.');
end

function data = set_deterministic_local_rhs_(data)
  % Ensure each sub(i).g exists, has correct length, and contains nontrivial
  % values at idx_c and idx_d.

  nSub = numel(data.sub);

  for i = 1:nSub
    if isfield(data.sub(i),'S') && ~isempty(data.sub(i).S)
      nGamma = size(data.sub(i).S, 1);
    else
      idx_all = [data.sub(i).idx_c(:); data.sub(i).idx_d(:)];
      nGamma = max([idx_all; 0]);
    end

    assert_true_(nGamma > 0, 'Could not determine local gamma size for subdomain %d.', i);

    gi = zeros(nGamma, 1);

    idx_d = data.sub(i).idx_d(:);
    if ~isempty(idx_d)
      assert_true_(max(idx_d) <= nGamma, 'sub(%d).idx_d out of bounds for g.', i);
      gi(idx_d) = (i) * (1:numel(idx_d))';
    end

    idx_c = data.sub(i).idx_c(:);
    if ~isempty(idx_c)
      assert_true_(max(idx_c) <= nGamma, 'sub(%d).idx_c out of bounds for g.', i);
      gi(idx_c) = (100 + 7*i) * ones(numel(idx_c), 1);
    end

    data.sub(i).g = gi;
  end
end

function [gC, gD] = assemble_rhs_pieces_like_solve_fetidp_(data)
  % Mirror the assembly in solve_fetidp.m.
  sub  = data.sub;
  nSub = numel(sub);

  gD = zeros(data.nDeltaProd, 1);
  gC = zeros(data.primal.nC, 1);

  assert_true_(isfield(data,'delta_range') && iscell(data.delta_range), 'Missing/invalid data.delta_range.');
  assert_true_(numel(data.delta_range) == nSub, 'data.delta_range must have one entry per subdomain.');

  for i = 1:nSub
    assert_true_(isfield(sub(i),'g'), 'sub(%d).g missing.', i);
    gi = sub(i).g(:);

    rng = data.delta_range{i};
    if ~isempty(rng)
      idx_d = sub(i).idx_d(:);
      assert_true_(numel(rng) == numel(idx_d), 'delta_range{%d} length must match numel(sub(%d).idx_d).', i, i);
      assert_true_(max(idx_d) <= numel(gi), 'sub(%d).idx_d out of bounds for sub(%d).g.', i, i);
      gD(rng) = gi(idx_d);
    end

    c_ids = sub(i).c_ids(:);
    if ~isempty(c_ids)
      idx_c = sub(i).idx_c(:);
      assert_true_(numel(c_ids) == numel(idx_c), 'sub(%d): numel(c_ids) must match numel(idx_c).', i);
      assert_true_(max(idx_c) <= numel(gi), 'sub(%d).idx_c out of bounds for sub(%d).g.', i, i);
      gC(c_ids) = gC(c_ids) + gi(idx_c);
    end
  end
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

  error('test_solve_fetidp:paths', 'could not locate project root containing main/setup_paths.m');
end

function assert_true_(cond, msg, varargin)
  if nargin < 2, msg = 'assertion failed'; end
  if ~cond
    if nargin > 2
      msg = sprintf(msg, varargin{:});
    end
    error('test_solve_fetidp:assert', '%s', msg);
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
    error('test_solve_fetidp:expected_error', 'Expected an error but none was thrown.');
  end
end