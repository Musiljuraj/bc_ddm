function test_applyM_lambda()
%TEST_APPLYM_LAMBDA  Unit tests for applyM_lambda().
%
% Location (suggested): tests/ch4_feti_dp/test_applyM_lambda.m
%
% This test prints a final PASS line if (and only if) all assertions pass.

  ensure_project_paths_();

  try
    % ---------------------------
    % Build a tiny synthetic case
    % ---------------------------
    nLambda = 3;
    nDelta  = 4;

    % Define Bd and BdT consistently.
    Bd  = [ 1, 0, 2, 0;
            0, 1, 0, 3;
           -1, 0, 0, 1 ];               % (nLambda x nDelta)
    BdT = Bd.';                          % (nDelta x nLambda)

    % Weights in packed-Delta space
    D = [1.0; 2.0; 0.5; 1.5];

    % Noncontiguous, disjoint ranges (includes empty range)
    delta_range = { [1, 3], [2, 4], [] };

    Sdd = cell(1,3);
    Sdd{1} = [2.0,  0.1;
             -0.3,  1.0];
    Sdd{2} = [1.0, -1.0;
              4.0,  0.5];
    Sdd{3} = []; % unused (empty range)

    data = struct();
    data.BdT          = BdT;
    data.Bd           = Bd;
    data.DeltaWeights = D;
    data.delta_range  = delta_range;
    data.Sdd          = Sdd;
    data.sub          = struct('id', {1,2,3}); % only used for numel(data.sub)

    % Sanity-check fixture invariants
    assert_true_(numel(data.DeltaWeights) == nDelta, 'DeltaWeights length mismatch in test fixture');
    assert_true_(size(data.BdT,1) == nDelta && size(data.BdT,2) == nLambda, 'BdT size mismatch in test fixture');
    assert_true_(size(data.Bd,1)  == nLambda && size(data.Bd,2)  == nDelta, 'Bd size mismatch in test fixture');
    assert_true_(ranges_are_disjoint_(data.delta_range, nDelta), 'delta_range must be disjoint and in-bounds for this test');

    % ---------------------------
    % Correctness vs explicit reference operator
    % ---------------------------
    r = [1.0; -2.0; 0.5];

    z = applyM_lambda(r, data);

    % Reference: z = Bd * (D .* (Sfull * (D .* (BdT*r))))
    t = data.BdT * r;
    t = data.DeltaWeights .* t;

    Sfull = zeros(nDelta, nDelta);
    for i = 1:numel(data.sub)
      rng = data.delta_range{i};
      if isempty(rng), continue; end
      Sfull(rng, rng) = data.Sdd{i};
    end

    uref = Sfull * t;
    uref = data.DeltaWeights .* uref;
    zref = data.Bd * uref;

    tol = 100 * eps(max(1, norm(zref)));
    assert_close_(z, zref, tol, 'applyM_lambda does not match explicit reference');

    % Output shape
    assert_true_(isvector(z) && rows(z) == nLambda && columns(z) == 1, 'z must be (nLambda x 1)');

    % ---------------------------
    % Linearity and shape handling
    % ---------------------------
    r1 = [1.0; 0.0; -1.0];
    r2 = [0.5; -1.5; 2.0];
    a  = -0.3;
    b  =  2.0;

    z1 = applyM_lambda(r1, data);
    z2 = applyM_lambda(r2, data);
    zlin = applyM_lambda(a*r1 + b*r2, data);

    assert_close_(zlin, a*z1 + b*z2, 100*eps(max(1, norm(a*z1 + b*z2))), 'linearity check failed');

    % Row-vector input should be accepted (routine reshapes with r(:))
    zrow = applyM_lambda(r1.', data);
    assert_close_(zrow, z1, 100*eps(max(1, norm(z1))), 'row-vector input handling failed');

    % ---------------------------
    % Invariance to subdomain ordering (disjoint ranges)
    % ---------------------------
    perm = [2, 1, 3];
    data2 = data;
    data2.sub         = data.sub(perm);
    data2.delta_range = data.delta_range(perm);
    data2.Sdd         = data.Sdd(perm);

    z_perm = applyM_lambda(r, data2);
    assert_close_(z_perm, z, 100*eps(max(1, norm(z))), 'subdomain reordering changed the result');

    % ---------------------------
    % Negative tests (must throw)
    % ---------------------------
    assert_throws_(@() applyM_lambda(), 'wrong arg count (0 args) did not throw');
    assert_throws_(@() applyM_lambda(r), 'wrong arg count (1 arg) did not throw');
    assert_throws_(@() applyM_lambda(r, data, 123), 'wrong arg count (3 args) did not throw');

    % Size mismatch: r length incompatible with BdT
    assert_throws_(@() applyM_lambda([1;2], data), 'size-mismatched r did not throw');

    % Missing required field
    data_missing = rmfield(data, 'Bd');
    assert_throws_(@() applyM_lambda(r, data_missing), 'missing data.Bd did not throw');

    % DeltaWeights length mismatch
    data_badD = data;
    data_badD.DeltaWeights = ones(nDelta-1, 1);
    assert_throws_(@() applyM_lambda(r, data_badD), 'DeltaWeights length mismatch did not throw');

    % Invalid residual types (should throw after hardening)
    bad_inputs = { struct('x',1), 'abc', [1+1i; 2; 3], [1; NaN; 3], [1; Inf; 3], true(3,1) };
    bad_tags   = { 'struct',       'char', 'complex',      'NaN',      'Inf',      'logical' };
    for k = 1:numel(bad_inputs)
      threw = did_throw_(@() applyM_lambda(bad_inputs{k}, data));
      if ~threw
        error('test_applyM_lambda:ExpectedErrorNotThrown', ...
              'Expected an error for bad input type (%s), but none was thrown.', bad_tags{k});
      end
    end

    % If we get here, everything passed.
    fprintf('test_applyM_lambda: PASS\n');

  catch err
    % Print a clear FAIL line, then rethrow to preserve failure semantics.
    fprintf('test_applyM_lambda: FAIL (%s)\n', err.message);
    rethrow(err);
  end

end

% ============================================================
% Local helpers
% ============================================================

function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did
    return;
  end

  % Strategy 1: setup_paths already resolvable
  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    main_dir = fileparts(sp);
    root_dir = fileparts(main_dir);
    addpath(genpath(fullfile(root_dir, 'tests')));
    did = true;
    return;
  end

  % Strategy 2: derive root from this test file location
  thisdir = fileparts(mfilename('fullpath'));
  root_dir = thisdir;

  found = false;
  for k = 1:10
    if exist(fullfile(root_dir, 'main', 'setup_paths.m'), 'file') == 2
      found = true;
      break;
    end
    parent = fileparts(root_dir);
    if isequal(parent, root_dir)
      break;
    end
    root_dir = parent;
  end

  if ~found
    error('ensure_project_paths_:rootNotFound', ...
          'could not locate project root containing main/setup_paths.m');
  end

  addpath(fullfile(root_dir, 'main'));
  setup_paths();
  addpath(genpath(fullfile(root_dir, 'tests')));

  did = true;
end

function assert_close_(a, b, tol, msg)
  if nargin < 4, msg = 'values differ'; end
  if any(size(a) ~= size(b))
    error('test_applyM_lambda:SizeMismatch', '%s (size mismatch)', msg);
  end
  if norm(a - b) > tol
    error('test_applyM_lambda:NotClose', '%s (||a-b||=%g, tol=%g)', msg, norm(a-b), tol);
  end
end

function assert_true_(cond, msg)
  if nargin < 2, msg = 'assertion failed'; end
  if ~cond
    error('test_applyM_lambda:AssertTrueFailed', '%s', msg);
  end
end

function assert_throws_(fh, msg)
  if nargin < 2, msg = 'expected an error but none was thrown'; end
  threw = did_throw_(fh);
  if ~threw
    error('test_applyM_lambda:ExpectedErrorNotThrown', '%s', msg);
  end
end

function threw = did_throw_(fh)
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
end

function ok = ranges_are_disjoint_(ranges, n)
  used = false(n,1);
  ok = true;
  for i = 1:numel(ranges)
    rng = ranges{i};
    if isempty(rng), continue; end
    if ~isvector(rng), ok = false; return; end
    rng = rng(:);
    if any(rng < 1) || any(rng > n), ok = false; return; end
    if any(abs(rng - round(rng)) > 0), ok = false; return; end
    if any(used(rng)), ok = false; return; end
    used(rng) = true;
  end
end