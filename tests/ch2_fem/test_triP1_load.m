function test_triP1_load()
%TEST_TRIP1_LOAD Strengthened unit tests for triP1_load.m (P1 element load vector).

  % ----------------------------
  % Path setup (robust to runner)
  % ----------------------------
  % CHANGED: only call setup_paths if it exists on the path
  if exist('setup_paths', 'file') == 2
    setup_paths();
  end

  % Small helper tolerance (scale-aware)
  % ADDED
  local_tol = @(scale) (50 * eps(max(1, scale)));

  % ----------------------------
  % Test data: reference triangle
  % ----------------------------
  xy_ref = [0 0; 1 0; 0 1];
  area_ref = 0.5;
  xc_ref = 1/3; yc_ref = 1/3;

  % ============================================================
  % Test 0: Wrong-arg-count must error
  % ============================================================
  % ADDED
  assert_throws(@() triP1_load(),                 'bad nargin case 0 did not throw');
  assert_throws(@() triP1_load(xy_ref),           'bad nargin case 1 did not throw');
  assert_throws(@() triP1_load(xy_ref, @(x,y) 1, 7), 'bad nargin case 3 did not throw');

  % ============================================================
  % Test 1: Constant f on reference triangle (exact for centroid rule)
  % ============================================================
  c = 2.5;
  fconst = @(x,y) c;

  fe = triP1_load(xy_ref, fconst);

  % ADDED: output shape/type invariants
  assert(isnumeric(fe), 'fe must be numeric');
  assert(size(fe,1) == 3 && size(fe,2) == 1, 'fe must be 3x1');

  fe_expected = (c * area_ref / 3) * ones(3,1);
  tol = local_tol(norm(fe_expected, inf));

  assert_allclose(fe, fe_expected, tol, 'Test 1 failed: constant f value mismatch');

  % Invariant: all entries equal for centroid rule
  % ADDED
  assert_allclose(fe(1), fe(2), local_tol(abs(fe(1))), 'Test 1 failed: entries not equal (1 vs 2)');
  assert_allclose(fe(2), fe(3), local_tol(abs(fe(2))), 'Test 1 failed: entries not equal (2 vs 3)');

  % Invariant: sum(fe) = f(xc,yc) * area
  assert_allclose(sum(fe), c * area_ref, local_tol(abs(c*area_ref)), ...
                 'Test 1 failed: sum(fe) does not match c*area');

  % ============================================================
  % Test 2: Permutation / orientation invariance (all 6 orderings)
  % ============================================================
  % CHANGED: strengthened from a single CW/CCW pair to all permutations
  perms_idx = perms([1 2 3]);  % 6x3
  fe0 = triP1_load(xy_ref, fconst);
  for k = 1:size(perms_idx,1)
    p = perms_idx(k,:);
    fe_p = triP1_load(xy_ref(p,:), fconst);
    tol = local_tol(norm(fe0, inf));
    assert_allclose(fe_p, fe0, tol, sprintf('Test 2 failed: permutation k=%d changed fe', k));
  end

  % ============================================================
  % Test 3: Centroid is used (non-constant f, compare to centroid formula)
  % ============================================================
  % ADDED
  f_lin = @(x,y) (x + 2*y + 0.1);  % simple deterministic function
  fval = f_lin(xc_ref, yc_ref);
  fe_expected = (fval * area_ref / 3) * ones(3,1);

  fe = triP1_load(xy_ref, f_lin);
  tol = local_tol(norm(fe_expected, inf));
  assert_allclose(fe, fe_expected, tol, 'Test 3 failed: centroid formula mismatch');

  % ============================================================
  % Test 4: Degenerate triangle must throw (collinear points)
  % ============================================================
  xy_deg = [0 0; 1 0; 2 0];
  assert_throws(@() triP1_load(xy_deg, fconst), 'Test 4 failed: degenerate triangle did not throw'); % FIXED msg

  % ============================================================
  % Test 5: Input shape validation for xy must throw
  % ============================================================
  % ADDED
  assert_throws(@() triP1_load([0 0; 1 0], fconst),          'Test 5 failed: 2x2 xy did not throw');
  assert_throws(@() triP1_load(zeros(3,3), fconst),          'Test 5 failed: 3x3 xy did not throw');
  assert_throws(@() triP1_load(zeros(3,1), fconst),          'Test 5 failed: 3x1 xy did not throw');

  % ============================================================
  % Test 6: f_handle must be a function handle
  % ============================================================
  % ADDED
  assert_throws(@() triP1_load(xy_ref, 3.14), 'Test 6 failed: non-handle f did not throw');

  % ============================================================
  % Test 7: f_handle return value must be finite scalar
  % ============================================================
  % ADDED
  assert_throws(@() triP1_load(xy_ref, @(x,y) [1;2]), 'Test 7 failed: vector fval did not throw');
  assert_throws(@() triP1_load(xy_ref, @(x,y) NaN),   'Test 7 failed: NaN fval did not throw');
  assert_throws(@() triP1_load(xy_ref, @(x,y) Inf),   'Test 7 failed: Inf fval did not throw');

  % ============================================================
  % Test 8: xy should be real, numeric, finite (contract hardening)
  % ============================================================
  % ADDED
  %
  % NOTE: These are expected to FAIL with the current triP1_load implementation,
  % because it does not validate xy for numeric/finite/real. Per project rules,
  % the preferred resolution is to harden triP1_load (Option 1), not to relax tests.
  %
  bad_xy = {
    [NaN 0; 1 0; 0 1],  'xy with NaN'
    [Inf 0; 1 0; 0 1],  'xy with Inf'
    complex(xy_ref),    'xy complex'
    char(xy_ref),       'xy char (ASCII coercion)'
    logical(xy_ref),    'xy logical (0/1 coercion)'
  };
  for i = 1:size(bad_xy,1)
    xy_bad = bad_xy{i,1};
    label  = bad_xy{i,2};
    assert_throws(@() triP1_load(xy_bad, fconst), ...
                  sprintf('Test 8 failed: expected error not thrown for case %d (%s)', i, label));
  end

  fprintf('PASS: test_triP1_load\n');
end


% ============================================================
% Helpers
% ============================================================

% ADDED
function assert_throws(fh, msg)
  if nargin < 2
    msg = 'Expected an error but none was thrown.';
  end
  did_error = false;
  try
    fh();
  catch
    did_error = true;
  end
  assert(did_error, msg);
end

% ADDED
function assert_allclose(a, b, tol, msg)
  if nargin < 4
    msg = 'Values are not close.';
  end
  err = norm(a - b, inf);
  assert(err <= tol, sprintf('%s (err=%g, tol=%g)', msg, err, tol));
end