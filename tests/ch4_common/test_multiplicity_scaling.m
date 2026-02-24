function test_multiplicity_scaling()
%TEST_MULTIPLICITY_SCALING  Unit tests for multiplicity_scaling(prod).
%
% Suggested location: tests/ch4_common/test_multiplicity_scaling.m
%
% Run from project root:
%   addpath('main'); setup_paths(); addpath(genpath('tests')); test_multiplicity_scaling();

  ensure_project_paths_();

  % -----------------------
  % Positive: simple groups
  % -----------------------
  prod = struct();
  prod.nProd = 4;
  prod.hat2prod = { [1 2], [3], [4] };

  omega = multiplicity_scaling(prod);

  assert(isnumeric(omega));
  assert(iscolumn(omega));
  assert(numel(omega) == prod.nProd);

  assert_close_(omega, [0.5; 0.5; 1.0; 1.0], 50*eps, 0, 'basic weighting mismatch');

  % Each hat group should sum to 1 (assuming disjoint partition semantics).
  for h = 1:numel(prod.hat2prod)
    idx = prod.hat2prod{h}(:);
    if isempty(idx), continue; end
    assert_close_(sum(omega(idx)), 1.0, 50*eps, 0, 'group sum must be 1');
  end

  % Permutation invariance within a hat group.
  prod2 = prod;
  prod2.hat2prod{1} = [2 1];
  omega2 = multiplicity_scaling(prod2);
  assert_close_(omega2, omega, 0, 0, 'permutation within group should not change omega');

  % -------------------------
  % Positive: another mapping
  % -------------------------
  prod = struct();
  prod.nProd = 6;
  prod.hat2prod = { [1 4 6], [2 3], [5] };

  omega = multiplicity_scaling(prod);

  expected = zeros(6,1);
  expected([1 4 6]) = 1/3;
  expected([2 3])   = 1/2;
  expected(5)       = 1;

  assert_close_(omega, expected, 50*eps, 0, 'second fixture weighting mismatch');

  % --------------------------------
  % Negative: wrong argument counts
  % --------------------------------
  assert_throws_(@() multiplicity_scaling(), 'nargin=0');
  assert_throws_(@() multiplicity_scaling(prod, 1), 'nargin=2');

  % -------------------------------
  % Negative: prod type / fields
  % -------------------------------
  assert_throws_(@() multiplicity_scaling(123), 'prod not struct');
  assert_throws_(@() multiplicity_scaling(struct()), 'missing fields');
  assert_throws_(@() multiplicity_scaling(struct('nProd',1)), 'missing hat2prod');
  assert_throws_(@() multiplicity_scaling(struct('hat2prod',{{[1]}})), 'missing nProd');

  % ------------------------
  % Negative: nProd invalid
  % ------------------------
  bad = prod; bad.nProd = -1;
  assert_throws_(@() multiplicity_scaling(bad), 'nProd negative');

  bad = prod; bad.nProd = 2.5;
  assert_throws_(@() multiplicity_scaling(bad), 'nProd nonint');

  bad = prod; bad.nProd = [2 3];
  assert_throws_(@() multiplicity_scaling(bad), 'nProd nonscalar');

  % ----------------------------
  % Negative: hat2prod invalid
  % ----------------------------
  bad = prod; bad.hat2prod = 42;
  assert_throws_(@() multiplicity_scaling(bad), 'hat2prod not cell');

  % --------------------------------------
  % Negative: uncovered product index -> error
  % --------------------------------------
  bad = struct('nProd', 3, 'hat2prod', {{ [1], [2] }});
  assert_throws_(@() multiplicity_scaling(bad), 'uncovered index');

  % --------------------------------------
  % Negative: out-of-range / non-integer indices
  % --------------------------------------
  bad = struct('nProd', 3, 'hat2prod', {{ [0], [1 2], [3] }});
  assert_throws_(@() multiplicity_scaling(bad), 'idx zero');

  bad = struct('nProd', 3, 'hat2prod', {{ [1], [2], [4] }});
  assert_throws_(@() multiplicity_scaling(bad), 'idx too large');

  bad = struct('nProd', 2, 'hat2prod', {{ [1.2], [2] }});
  assert_throws_(@() multiplicity_scaling(bad), 'idx nonint');

  % --------------------------------------
  % Negative: duplicates/overlaps should error
  % (Strong partition validity requirements.)
  % --------------------------------------
  bad = struct('nProd', 2, 'hat2prod', {{ [1 1], [2] }});
  assert_throws_(@() multiplicity_scaling(bad), 'dup within group');

  bad = struct('nProd', 3, 'hat2prod', {{ [1 2], [2 3] }});
  assert_throws_(@() multiplicity_scaling(bad), 'overlap across groups');

  fprintf('PASS: test_multiplicity_scaling\n');
end

% =========================
% Helpers (mirrors test_apply_local_schur.m path logic)
% =========================

function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did
    return;
  end

  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    if ~isempty(sp)
      maindir = fileparts(sp);
      rootdir = fileparts(maindir);
      addpath(genpath(fullfile(rootdir, 'tests')));
    end
    did = true;
    return;
  end

  thisdir  = fileparts(mfilename('fullpath'));
  testsdir = fileparts(thisdir);
  rootdir  = fileparts(testsdir);
  maindir  = fullfile(rootdir, 'main');

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, 'tests')));

  did = true;
end

function assert_close_(a, b, reltol, abstol, msg)
  if nargin < 3 || isempty(reltol), reltol = 1e-12; end
  if nargin < 4 || isempty(abstol), abstol = 0; end
  if nargin < 5, msg = 'values not close'; end

  assert(isequal(size(a), size(b)));

  da = double(a(:));
  db = double(b(:));

  scale = max(1.0, norm(db, inf));
  err = norm(da - db, inf);

  if err > max(abstol, reltol * scale)
    error('test_multiplicity_scaling:NotClose', '%s (err=%g, scale=%g)', msg, err, scale);
  end
end

function assert_throws_(fh, casename)
  if nargin < 2, casename = ''; end
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  if ~threw
    if ~isempty(casename)
      fprintf(2, 'Expected an error but none was thrown (case: %s)\n', casename);
    end
    error('test_multiplicity_scaling:ExpectedError', 'Expected an error but none was thrown.');
  end
end