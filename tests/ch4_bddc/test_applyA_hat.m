function test_applyA_hat()
%TEST_APPLYA_HAT  Unit tests for applyA_hat (path-robust + non-interactive cleanup).
%
% This test prints an explicit PASS line on successful completion.

  fprintf("test_applyA_hat ... ");

  ensure_project_paths_();

  orig_path = path();
  mock_dir  = tempname();
  mkdir(mock_dir);

  % Install mock apply_blockdiag_S at highest precedence.
  write_mock_apply_blockdiag_S_(mock_dir);
  addpath(mock_dir, "-begin");
  clear functions;  % ensure applyA_hat sees the mock

  have_confirm = (exist("confirm_recursive_rmdir", "builtin") == 5) || ...
                 (exist("confirm_recursive_rmdir", "file") == 2);

  c = onCleanup(@() cleanup_(orig_path, mock_dir, have_confirm)); %#ok<NASGU>

  run_positive_cases_();
  run_negative_cases_();

  fprintf("PASS\n");
end


% ------------------------- Positive cases -------------------------

function run_positive_cases_()
  global BLOCKDIAG_S_LAST_SUB BLOCKDIAG_S_LAST_W BLOCKDIAG_S_MATRIX;

  % Tiny deterministic operator sizes
  nHat  = 3;
  nProd = 4;

  R = [ 1  0  2;
        0  1 -1;
        1  1  0;
       -1  0  1 ];          % (nProd x nHat)

  % Symmetric product-space operator for sanity checks
  S = [ 4  1  0  0;
        1  3  0  0;
        0  0  2  0;
        0  0  0  5 ];       % (nProd x nProd)

  sub = struct("tag", 7);   % carried through to mock
  data = struct("sub", sub, "bddc", struct("R", R));

  % Case 1: basic correctness with column x
  x = [ 2; -1; 3 ];
  BLOCKDIAG_S_MATRIX   = S;
  BLOCKDIAG_S_LAST_SUB = [];
  BLOCKDIAG_S_LAST_W   = [];

  y = applyA_hat(x, data);

  expected = R' * (S * (R * x(:)));
  assert_close_(y, expected, 0, 0);

  % Output is column of correct length and finite
  assert(isequal(size(y), [nHat, 1]));
  assert(all(isfinite(y(:))));

  % Verify the mock saw exactly what applyA_hat should pass
  assert(isstruct(BLOCKDIAG_S_LAST_SUB));
  assert(BLOCKDIAG_S_LAST_SUB.tag == sub.tag);
  assert_close_(BLOCKDIAG_S_LAST_W, R * x(:), 0, 0);

  % Case 2: row-vector input should be treated equivalently (x(:) behavior)
  xr = [ 2, -1, 3 ];
  yr = applyA_hat(xr, data);
  assert_close_(yr, expected, 0, 0);
  assert(isequal(size(yr), [nHat, 1]));

  % Case 3: linearity check (invariant property)
  x1 = [ 1; 4; -2 ];
  x2 = [ -3; 0;  5 ];
  a  = -2.5;
  ylin1 = applyA_hat(a*x1 + x2, data);
  ylin2 = a*applyA_hat(x1, data) + applyA_hat(x2, data);
  assert_close_(ylin1, ylin2, 1e-12, 1e-14);
end


% ------------------------- Negative cases -------------------------

function run_negative_cases_()
  % Minimal consistent data
  R    = [1 0; 0 1; 1 1];     % (3 x 2)
  data = struct("sub", struct(), "bddc", struct("R", R));

  % 1) Wrong argument count must throw
  must_throw_named_(1, "wrong arg count: missing data", @() applyA_hat(1));
  must_throw_named_(2, "wrong arg count: extra arg",    @() applyA_hat(1, data, 7));

  % 2) data missing bddc.R must throw (explicit check in applyA_hat)
  must_throw_named_(3, "data missing bddc.R", @() applyA_hat([1;2], struct("sub", struct(), "bddc", struct())));

  % 3) Missing sub should throw (field access before mock call)
  must_throw_named_(4, "data missing sub", @() applyA_hat([1;2], struct("bddc", struct("R", R))));

  % 4) Dimension mismatch must throw (R*x)
  must_throw_named_(5, "dimension mismatch R*x", @() applyA_hat([1;2;3], data));
end


% ------------------------- Helpers -------------------------

function must_throw_named_(idx, name, fh)
  didThrow = false;
  try
    fh();
  catch
    didThrow = true;
  end
  if ~didThrow
    error("test_applyA_hat:ExpectedErrorNotThrown", ...
          "Expected an error but none was thrown for bad case %d (%s).", idx, name);
  end
end

function assert_close_(a, b, rtol, atol)
  if nargin < 3, rtol = 1e-12; end
  if nargin < 4, atol = 1e-14; end
  da = a - b;
  na = norm(da(:), 2);
  nb = norm(b(:), 2);
  scale = max(atol, rtol * max(1, nb));
  assert(na <= scale);
end

function write_mock_apply_blockdiag_S_(mock_dir)
  % Creates a mock apply_blockdiag_S.m in a temporary directory so applyA_hat
  % becomes deterministic in this unit test.
  fn = fullfile(mock_dir, "apply_blockdiag_S.m");
  fid = fopen(fn, "w");
  if fid < 0
    error("test_applyA_hat:IO", "Could not create mock apply_blockdiag_S.m");
  end
  fprintf(fid, "%s\n", "function Sw = apply_blockdiag_S(sub, w)");
  fprintf(fid, "%s\n", "  global BLOCKDIAG_S_LAST_SUB BLOCKDIAG_S_LAST_W BLOCKDIAG_S_MATRIX;");
  fprintf(fid, "%s\n", "  BLOCKDIAG_S_LAST_SUB = sub;");
  fprintf(fid, "%s\n", "  BLOCKDIAG_S_LAST_W   = w;");
  fprintf(fid, "%s\n", "  S = BLOCKDIAG_S_MATRIX;");
  fprintf(fid, "%s\n", "  Sw = S * w;");
  fprintf(fid, "%s\n", "end");
  fclose(fid);
end

function cleanup_(orig_path, mock_dir, have_confirm)
  % Restore path and clear caches first.
  path(orig_path);
  clear apply_blockdiag_S;
  clear functions;

  % Make recursive removal non-interactive in Octave (if available).
  if have_confirm
    old_confirm = confirm_recursive_rmdir(false, "local");
  else
    old_confirm = [];
  end

  if exist(mock_dir, "dir") == 7
    rmdir(mock_dir, "s");
  end

  if have_confirm
    confirm_recursive_rmdir(old_confirm, "local");
  end
end

function ensure_project_paths_()
  % Required ensure-path mechanism (idempotent).
  persistent did;
  if ~isempty(did) && did
    return;
  end

  % Strategy 1: setup_paths already resolvable
  if exist("setup_paths", "file") == 2
    setup_paths();
    sp = which("setup_paths");
    main_dir = fileparts(sp);
    root_dir = fileparts(main_dir);
    addpath(genpath(fullfile(root_dir, "tests")));
    did = true;
    return;
  end

  % Strategy 2: derive root from this test file location
  thisdir = fileparts(mfilename("fullpath"));
  root_dir = thisdir;

  for k = 1:10
    candidate = fullfile(root_dir, "main", "setup_paths.m");
    if exist(candidate, "file") == 2
      addpath(fullfile(root_dir, "main"));
      setup_paths();
      addpath(genpath(fullfile(root_dir, "tests")));
      did = true;
      return;
    end
    parent = fileparts(root_dir);
    if strcmp(parent, root_dir)
      break;
    end
    root_dir = parent;
  end

  error("test_applyA_hat:RootNotFound", ...
        "could not locate project root containing main/setup_paths.m");
end