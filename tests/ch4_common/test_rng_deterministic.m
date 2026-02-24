function test_rng_deterministic()
%TEST_RNG_DETERMINISTIC  Unit tests for rng_deterministic(seed).
%
% Covers:
%   - Determinism for rand and randn for a fixed seed
%   - Different seeds produce different sequences (deterministic check)
%   - Signature checks (nargin/nargout)
%   - Invalid input checks (shape/type/value)

  ensure_project_paths_();

  % Preserve current RNG state (best-effort) so we don't pollute other tests.
  old_seed_rand  = [];
  old_seed_randn = [];
  try old_seed_rand  = rand('seed');  catch, end
  try old_seed_randn = randn('seed'); catch, end

  unwind_protect

    % ----------------------------
    % Case A: same seed => identical sequences (rand and randn)
    % ----------------------------
    rng_deterministic(7);
    a1 = rand(1, 8);
    b1 = randn(1, 8);
    a2 = rand(1, 8);
    b2 = randn(1, 8);

    rng_deterministic(7);
    a1b = rand(1, 8);
    b1b = randn(1, 8);
    a2b = rand(1, 8);
    b2b = randn(1, 8);

    assert(isequal(a1, a1b), 'rand sequence not reproducible for same seed.');
    assert(isequal(b1, b1b), 'randn sequence not reproducible for same seed.');
    assert(isequal(a2, a2b), 'rand sequence (continued) not reproducible for same seed.');
    assert(isequal(b2, b2b), 'randn sequence (continued) not reproducible for same seed.');

    % ----------------------------
    % Case B: different seeds => different sequences (deterministic check)
    % ----------------------------
    rng_deterministic(7);
    v7r = rand(1, 12);
    v7n = randn(1, 12);

    rng_deterministic(8);
    v8r = rand(1, 12);
    v8n = randn(1, 12);

    assert(~isequal(v7r, v8r) || ~isequal(v7n, v8n), ...
           'Different seeds unexpectedly produced identical sequences.');

    % ----------------------------
    % Signature checks
    % ----------------------------
    assert_throws_(@() rng_deterministic(),     'rng_deterministic nargin=0');
    assert_throws_(@() rng_deterministic(1,2),  'rng_deterministic nargin=2');

    % Requesting an output must throw (function defines no outputs)
    assert_throws_(@() eval('x = rng_deterministic(1);'), 'rng_deterministic nargout requested'); %#ok<NASGU>

    % ----------------------------
    % Invalid inputs must throw
    % ----------------------------
    assert_throws_(@() rng_deterministic([1 2]),   'non-scalar seed');
    assert_throws_(@() rng_deterministic(NaN),     'NaN seed');
    assert_throws_(@() rng_deterministic(Inf),     'Inf seed');
    assert_throws_(@() rng_deterministic(-Inf),    '-Inf seed');
    assert_throws_(@() rng_deterministic(1.25),    'non-integer seed');

    % Complex should not be accepted (exact failure mode may vary, but must throw)
    assert_throws_(@() rng_deterministic(1 + 2i),  'complex seed');

    % Per hardening baseline: reject logical and char (if these pass, it’s a real weakness)
    assert_throws_(@() rng_deterministic(true),    'logical seed');
    assert_throws_(@() rng_deterministic('7'),     'char seed');

    fprintf('PASS: test_rng_deterministic\n');

  unwind_protect_cleanup
    % Restore RNG state (best-effort)
    try
      if ~isempty(old_seed_rand)
        rand('seed', old_seed_rand);
      end
    catch
    end
    try
      if ~isempty(old_seed_randn)
        randn('seed', old_seed_randn);
      end
    catch
    end
  end_unwind_protect

end

% =========================
% Helpers
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

function assert_throws_(fh, label)
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  if ~threw
    fprintf('DIAG: Expected an error but none was thrown for case: %s\n', label);
  end
  assert(threw, sprintf('Expected an error, but none was thrown (%s).', label));
end