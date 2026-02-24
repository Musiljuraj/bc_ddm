function test_applyA_lambda()
%TEST_APPLYA_LAMBDA  Unit tests for applyA_lambda (path-robust + non-interactive cleanup).
%
% This test prints an explicit PASS line on successful completion.

  fprintf("test_applyA_lambda ... ");

  ensure_project_paths_();

  orig_path = path();
  mock_dir = tempname();
  mkdir(mock_dir);

  % Install mock solve_tildeS at highest precedence.
  write_mock_solve_tildeS_(mock_dir);
  addpath(mock_dir, "-begin");
  clear functions;  % ensure applyA_lambda sees the mock

  old_confirm = [];
  have_confirm = (exist("confirm_recursive_rmdir", "builtin") == 5) || ...
                 (exist("confirm_recursive_rmdir", "file") == 2);

  unwind_protect
    run_positive_cases_();
    run_negative_cases_();

    % If we reach here, everything passed.
    fprintf("PASS\n");  % ADDED
  unwind_protect_cleanup
    % Restore path and clear caches first.
    path(orig_path);
    clear functions;

    % Make recursive removal non-interactive in Octave (if available).
    if have_confirm
      old_confirm = confirm_recursive_rmdir(false, "local");
    end

    if exist(mock_dir, "dir") == 7
      rmdir(mock_dir, "s");
    end

    if have_confirm
      confirm_recursive_rmdir(old_confirm, "local");
    end
  end_unwind_protect
end


% ------------------------- Positive cases -------------------------

function run_positive_cases_()
  global SOLVE_TILDES_LAST_RC SOLVE_TILDES_LAST_T SOLVE_TILDES_LAST_DATA;

  % Tiny deterministic operator sizes
  nLambda = 2;
  nDelta  = 3;

  BdT = [ 1  2;
          0  1;
         -1  0 ];         % (nDelta x nLambda)

  Bd  = [ 1  0  1;
          0  2 -1 ];      % (nLambda x nDelta)

  data = struct("BdT", BdT, "Bd", Bd);

  % Case 1: basic correctness with column x
  x = [ 3; -2 ];
  SOLVE_TILDES_LAST_RC = [];
  SOLVE_TILDES_LAST_T  = [];
  SOLVE_TILDES_LAST_DATA = [];

  y = applyA_lambda(x, data);

  expected = Bd * (BdT * x(:));  % mock solve_tildeS returns wD = t
  assert_close_(y, expected, 0, 0);

  % Output is column of correct length
  assert(isequal(size(y), [nLambda, 1]));

  % Verify the mock saw exactly what applyA_lambda should pass
  assert(isempty(SOLVE_TILDES_LAST_RC));  % applyA_lambda currently passes []
  assert_close_(SOLVE_TILDES_LAST_T, BdT * x(:), 0, 0);
  assert(isstruct(SOLVE_TILDES_LAST_DATA));
  assert_close_(SOLVE_TILDES_LAST_DATA.BdT, BdT, 0, 0);
  assert_close_(SOLVE_TILDES_LAST_DATA.Bd,  Bd,  0, 0);

  % Case 2: row-vector input should be treated equivalently (x(:) behavior)
  xr = [ 3, -2 ];
  yr = applyA_lambda(xr, data);
  assert_close_(yr, expected, 0, 0);
  assert(isequal(size(yr), [nLambda, 1]));

  % Case 3: linearity check (invariant property)
  x1 = [ 1; 4 ];
  x2 = [ -2; 5 ];
  a  = -3;
  ylin1 = applyA_lambda(a*x1 + x2, data);
  ylin2 = a*applyA_lambda(x1, data) + applyA_lambda(x2, data);
  assert_close_(ylin1, ylin2, 1e-12, 1e-14);
end


% ------------------------- Negative cases -------------------------

function run_negative_cases_()
  % Minimal consistent data for argument-count and type checks
  BdT = [ 1 0;
          0 1;
          1 1 ];
  Bd  = [ 1 0 0;
          0 1 1 ];
  data = struct("BdT", BdT, "Bd", Bd);

  % 1) Wrong argument count must throw
  must_throw_named_(1, "wrong arg count: missing data", @() applyA_lambda(1));
  must_throw_named_(2, "wrong arg count: extra arg",    @() applyA_lambda(1, data, 7));

  % 2) data missing fields must throw
  must_throw_named_(3, "data missing BdT", @() applyA_lambda([1;2], struct("Bd", Bd)));
  must_throw_named_(4, "data missing Bd",  @() applyA_lambda([1;2], struct("BdT", BdT)));

  % 3) Dimension mismatch must throw
  bad = struct("BdT", BdT(:,1), "Bd", Bd);   % BdT wrong width
  must_throw_named_(5, "BdT dimension mismatch", @() applyA_lambda([1;2], bad));

  % 4) Strengthened input-contract expectations
  must_throw_named_(6, "x must be a vector (matrix rejected)", @() applyA_lambda(eye(2), data));
  must_throw_named_(7, "x must be numeric (char rejected)",    @() applyA_lambda("ab", data));
  must_throw_named_(8, "x must be real (complex rejected)",    @() applyA_lambda([1+1i; 2], data));
  must_throw_named_(9, "x must be finite (NaN rejected)",      @() applyA_lambda([NaN; 1], data));
  must_throw_named_(10,"x must be finite (Inf rejected)",      @() applyA_lambda([Inf; 1], data));
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
    error("test_applyA_lambda:ExpectedErrorNotThrown", ...
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

function write_mock_solve_tildeS_(mock_dir)
  % Creates a mock solve_tildeS.m in a temporary directory so applyA_lambda
  % becomes deterministic in this unit test.
  fn = fullfile(mock_dir, "solve_tildeS.m");
  fid = fopen(fn, "w");
  if fid < 0
    error("test_applyA_lambda:IO", "Could not create mock solve_tildeS.m");
  end
  fprintf(fid, "%s\n", "function out = solve_tildeS(rc, t, data)");
  fprintf(fid, "%s\n", "  global SOLVE_TILDES_LAST_RC SOLVE_TILDES_LAST_T SOLVE_TILDES_LAST_DATA;");
  fprintf(fid, "%s\n", "  SOLVE_TILDES_LAST_RC = rc;");
  fprintf(fid, "%s\n", "  SOLVE_TILDES_LAST_T  = t;");
  fprintf(fid, "%s\n", "  SOLVE_TILDES_LAST_DATA = data;");
  fprintf(fid, "%s\n", "  out = struct();");
  fprintf(fid, "%s\n", "  out.wD = t;  % identity map for deterministic testing");
  fprintf(fid, "%s\n", "end");
  fclose(fid);
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

  error("test_applyA_lambda:RootNotFound", ...
        "could not locate project root containing main/setup_paths.m");
end