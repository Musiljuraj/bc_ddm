% ============================================================
% File: tests/ch6_spectra/test_assemble_from_apply.m
% ============================================================
function test_assemble_from_apply()
%TEST_ASSEMBLE_FROM_APPLY_IDENTITY  Unit test for assemble_from_apply on identity operator.
%
% This test follows the same conventions as the existing test suite:
%   - it is path-robust (works regardless of current working directory),
%   - it prints an explicit PASS line on success,
%   - it uses small deterministic checks and clear error messages.
%
% Why "identity" in the name?
%   We test the most basic and unambiguous scenario: T(x) = x.
%   Then assembling by unit vectors e_j must return the identity matrix:
%     Mat(:,j) = T(e_j) = e_j  =>  Mat = I.
%   The scenario label in the filename makes the intent obvious and allows
%   adding other scenario tests later (e.g. diagonal map, error handling),
%   without overloading a single generic test file.

  fprintf("test_assemble_from_apply_identity ... ");

  ensure_project_paths_();

  % ----------------------------
  % Case 1: n = 1 (smallest nontrivial)
  % ----------------------------
  n = 1;
  applyFun = @(x) x;  % identity mapping
  Mat = assemble_from_apply(applyFun, n);
  assert_close_(Mat, eye(n), 0, 0, "Identity assembly failed for n=1.");

  % ----------------------------
  % Case 2: moderate n
  % ----------------------------
  n = 5;
  Mat = assemble_from_apply(applyFun, n);
  assert_close_(Mat, eye(n), 0, 0, "Identity assembly failed for n=5.");

  % ----------------------------
  % Case 3: another size + opts path (guards against off-by-one logic)
  % ----------------------------
  n = 17;
  opts = struct("verbose", false, "progress_every", 5, "force_full", true);
  Mat = assemble_from_apply(applyFun, n, opts);
  assert_close_(Mat, eye(n), 0, 0, "Identity assembly failed for n=17.");

  % ----------------------------
  % Case 4: wrong output length must error (dimension mismatch guard)
  % ----------------------------
  n = 4;
  badApply = @(x) [x; 0];  % returns length n+1
  assert_throws_(@() assemble_from_apply(badApply, n), ...
                 "InvalidApplyOutputSize");

  % If we reach here, everything passed.
  fprintf("PASS\n");
end


% ============================================================
% Helpers (consistent with existing suite patterns)
% ============================================================

function ensure_project_paths_()
% Ensure src/ and tests/ are on path regardless of where the test is run from.
% This matches the robust pattern used in ch3/ch4 tests.
  persistent did;
  if ~isempty(did) && did
    return;
  end

  if exist("setup_paths", "file") == 2
    setup_paths();
    sp = which("setup_paths");
    if ~isempty(sp)
      maindir = fileparts(sp);
      rootdir = fileparts(maindir);
      addpath(genpath(fullfile(rootdir, "tests")));
    end
    did = true;
    return;
  end

  % Fallback: locate project root from this test file location
  thisdir  = fileparts(mfilename("fullpath")); % .../tests/ch6_spectra
  testsdir = fileparts(thisdir);               % .../tests
  rootdir  = fileparts(testsdir);              % project root
  maindir  = fullfile(rootdir, "main");        % .../main

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, "tests")));

  did = true;
end

function assert_close_(a, b, atol, rtol, msg)
%ASSERT_CLOSE_  Assert a and b are close with absolute/relative tolerances.
  if nargin < 5, msg = "Values not close."; end
  if nargin < 4 || isempty(rtol), rtol = 0; end
  if nargin < 3 || isempty(atol), atol = 0; end

  if ~isequal(size(a), size(b))
    error("%s Size mismatch: got %dx%d, expected %dx%d.", msg, size(a,1), size(a,2), size(b,1), size(b,2));
  end

  da = a(:) - b(:);
  err = max(abs(da));
  tol = atol + rtol * max(1, max(abs(b(:))));
  assert(err <= tol, sprintf("%s Max abs error = %.3e (tol = %.3e).", msg, err, tol));
end

function assert_throws_(fhandle, must_contain)
%ASSERT_THROWS_  Assert fhandle throws; optionally require substring in message.
  did_throw = false;
  try
    fhandle();
  catch err
    did_throw = true;
    if nargin >= 2 && ~isempty(must_contain)
      assert(~isempty(strfind(err.message, must_contain)), ...
             "Exception message did not contain expected substring.");
    end
  end
  assert(did_throw, "Expected function to throw, but it did not.");
end