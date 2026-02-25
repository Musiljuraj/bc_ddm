function test_applyM_bddc()
%TEST_APPLYM_BDDC  Unit tests for applyM_bddc().
%
% Location (suggested): tests/ch4_bddc/test_applyM_bddc.m
%
% This test prints a final PASS line if (and only if) all assertions pass.

  ensure_project_paths_();

  try
    % ============================================================
    % Case 1: Identity-like fixture (no coarse, single subdomain)
    % Expect: z_hat == r_hat
    % ============================================================
    nHat  = 3;
    nProd = 3;

    data1 = struct();
    data1.sub = struct('prod_idx_d', { (1:nProd) });

    bddc1 = struct();
    bddc1.R     = eye(nProd, nHat);
    bddc1.omega = ones(nProd, 1);
    bddc1.Psi   = [];     % unused when K0_R is empty
    bddc1.K0_R  = [];     % no coarse correction
    Sdd1        = eye(nProd);
    bddc1.Sdd_R = { chol(Sdd1) };

    data1.bddc = bddc1;

    r1 = [1.0; -2.0; 0.5];

    z1 = applyM_bddc(r1, data1);
    tol = 100 * eps(max(1, norm(r1)));
    assert_close_(z1, r1, tol, 'identity-like fixture failed (z_hat should equal r_hat)');

    % Output shape + finiteness
    assert_true_(isvector(z1) && size(z1,1) == nHat && size(z1,2) == 1, 'z_hat must be (nHat x 1)');
    assert_true_(all(isfinite(z1)), 'z_hat must not contain NaN/Inf for valid inputs');

    % Row-vector input should be accepted (routine reshapes with r(:))
    z1row = applyM_bddc(r1.', data1);
    assert_close_(z1row, z1, 100*eps(max(1, norm(z1))), 'row-vector input handling failed');

    % ============================================================
    % Case 2: Nontrivial fixture with coarse correction (two subs)
    % Compare to a tiny explicit reference expression.
    % ============================================================
    data2 = struct();
    data2.sub = struct('prod_idx_d', { [1;2], 3 });

    bddc2 = struct();
    bddc2.R     = eye(nProd, nHat);
    bddc2.omega = [1; 2; 3];

    % Local solves on disjoint index sets
    SddA = diag([4, 5]);
    SddB = 6;
    bddc2.Sdd_R = { chol(SddA), chol(SddB) };

    % Coarse (rank-1) component: Psi * inv(Psi'*Psi) * Psi' * r_prod
    bddc2.Psi  = ones(nProd, 1);
    K0         = (bddc2.Psi' * bddc2.Psi); % = 3
    bddc2.K0_R = chol(K0);

    data2.bddc = bddc2;

    r = [1.0; -2.0; 0.5];
    z = applyM_bddc(r, data2);

    % Explicit reference:
    % r_prod = omega .* (R*r)
    r_prod = bddc2.omega .* (bddc2.R * r);

    z_sub = zeros(nProd,1);
    % sub 1: idx [1,2]
    idx = data2.sub(1).prod_idx_d(:);
    Rdd = bddc2.Sdd_R{1};
    z_sub(idx) = Rdd \ (Rdd' \ r_prod(idx));
    % sub 2: idx [3]
    idx = data2.sub(2).prod_idx_d(:);
    Rdd = bddc2.Sdd_R{2};
    z_sub(idx) = Rdd \ (Rdd' \ r_prod(idx));

    rhs0  = full(bddc2.Psi' * r_prod);
    alpha = bddc2.K0_R \ (bddc2.K0_R' \ rhs0);
    z0    = bddc2.Psi * alpha;

    z_prod_ref = z_sub + z0;
    z_ref = bddc2.R' * (bddc2.omega .* z_prod_ref);

    tol2 = 200 * eps(max(1, norm(z_ref)));
    assert_close_(z, z_ref, tol2, 'nontrivial fixture does not match explicit reference');

    % ============================================================
    % Negative tests (must throw)
    % ============================================================
    assert_throws_(@() applyM_bddc(), 'wrong arg count (0 args) did not throw');
    assert_throws_(@() applyM_bddc(r), 'wrong arg count (1 arg) did not throw');
    assert_throws_(@() applyM_bddc(r, data2, 123), 'wrong arg count (3 args) did not throw');

    data_missing = rmfield(data2, 'bddc');
    assert_throws_(@() applyM_bddc(r, data_missing), 'missing data.bddc did not throw');

    % Size mismatch: r length incompatible with R
    assert_throws_(@() applyM_bddc([1;2], data2), 'size-mismatched r_hat did not throw');

    % If we get here, everything passed.
    fprintf('test_applyM_bddc: PASS\n');

  catch err
    % Print a clear FAIL line, then rethrow to preserve failure semantics.
    fprintf('test_applyM_bddc: FAIL (%s)\n', err.message);
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
    error('test_applyM_bddc:SizeMismatch', '%s (size mismatch)', msg);
  end
  if norm(a - b) > tol
    error('test_applyM_bddc:NotClose', '%s (||a-b||=%g, tol=%g)', msg, norm(a-b), tol);
  end
end

function assert_true_(cond, msg)
  if nargin < 2, msg = 'assertion failed'; end
  if ~cond
    error('test_applyM_bddc:AssertTrueFailed', '%s', msg);
  end
end

function assert_throws_(fh, msg)
  if nargin < 2, msg = 'expected an error but none was thrown'; end
  threw = did_throw_(fh);
  if ~threw
    error('test_applyM_bddc:ExpectedErrorNotThrown', '%s', msg);
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