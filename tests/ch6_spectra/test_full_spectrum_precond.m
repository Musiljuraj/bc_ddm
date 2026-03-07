% ============================================================
% File: tests/ch6_spectra/test_full_spectrum_precond.m
% ============================================================
function test_full_spectrum_precond()
%TEST_FULL_SPECTRUM_PRECOND  Unit test for full_spectrum_precond (synthetic SPD case).
%
% Conventions (consistent with the Chapter 6 tests):
%   - Robust path setup (works regardless of current working directory)
%   - Concise PASS/FAIL communication
%   - Deterministic random seed
%   - No dependency on PDE/problem setup
%
% Test idea:
%   Generate small SPD matrices A and Minv (representing M^{-1}). Define
%   matrix-free actions applyA(x)=A*x and applyMinv(x)=Minv*x, call
%   full_spectrum_precond, and compare returned eigenvalues with eig(Minv*A).
%
% Correctness basis:
%   full_spectrum_precond forms K = R*A*R' with chol(Minv)=R (R'*R=Minv).
%   K is similar to Minv*A, hence they have identical spectra (in exact arithmetic).

  fprintf('test_full_spectrum_precond ... ');

  ensure_project_paths_();

  % ----------------------------
  % Deterministic RNG (reproducible test)
  % ----------------------------
  rand('seed', 1);
  randn('seed', 1);

  % ----------------------------
  % Build synthetic SPD matrices
  % ----------------------------
  n = 12;

  % SPD operator A = Q'*Q + alpha*I
  Q = randn(n, n);
  A = Q' * Q + 0.5 * eye(n);

  % SPD Minv = P'*P + beta*I
  P = randn(n, n);
  Minv = P' * P + 0.7 * eye(n);

  % Matrix-free apply routines (match the signature used in Chapter 6)
  applyA    = @(x) A * x;
  applyMinv = @(x) Minv * x;

  % Options: keep defaults but be explicit (readability)
  opts = struct();
  opts.symmetrize     = true;
  opts.store_matrices = false;
  opts.verbose        = false;

  % ----------------------------
  % Run function under test
  % ----------------------------
  spec = full_spectrum_precond(applyA, applyMinv, n, opts);

  % ----------------------------
  % Basic checks
  % ----------------------------
  assert(spec.chol_ok, ...
         'Expected chol(Minv) to succeed for synthetic SPD Minv.');

  assert(isfield(spec, 'eigvals'), ...
         'Expected spec.eigvals field to exist.');

  assert(numel(spec.eigvals) == n, ...
         'Expected spec.eigvals to contain n eigenvalues.');

  % ----------------------------
  % Reference eigenvalues: eig(Minv*A)
  % ----------------------------
  ref = eig(Minv * A);
  ref = real(ref(:));
  ref = sort(ref, 'ascend');

  got = spec.eigvals(:);

  % ----------------------------
  % Compare with tolerance
  % ----------------------------
  % Note: eig can be sensitive; use a modest but strict tolerance.
  assert_close_vec_(got, ref, 1e-9, 1e-9, ...
                    'Eigenvalues do not match reference eig(Minv*A).');

  fprintf('PASS\n');
end


% ============================================================
% Helpers (path setup + assertions)
% ============================================================

function ensure_project_paths_()
%ENSURE_PROJECT_PATHS_  Ensure src/ and tests/ are on path.
%
% This helper mirrors the robust approach used in the Chapter 6 tests:
% - If setup_paths is already visible, call it and also add tests/.
% - Otherwise, infer the project root from this test file location, add main/,
%   call setup_paths, and add tests/.

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

  % Fallback: locate project root from this test file location
  thisdir  = fileparts(mfilename('fullpath')); % .../tests/ch6_spectra
  testsdir = fileparts(thisdir);               % .../tests
  rootdir  = fileparts(testsdir);              % project root
  maindir  = fullfile(rootdir, 'main');        % .../main

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, 'tests')));

  did = true;
end


function assert_close_vec_(a, b, atol, rtol, msg)
%ASSERT_CLOSE_VEC_  Assert two vectors are close with abs/rel tolerance.

  if nargin < 5
    msg = 'Vectors not close.';
  end
  if nargin < 4 || isempty(rtol)
    rtol = 0;
  end
  if nargin < 3 || isempty(atol)
    atol = 0;
  end

  a = a(:);
  b = b(:);

  if numel(a) ~= numel(b)
    error('%s Length mismatch: got %d, expected %d.', ...
          msg, numel(a), numel(b));
  end

  err = max(abs(a - b));
  tol = atol + rtol * max(1, max(abs(b)));

  assert(err <= tol, ...
         sprintf('%s Max abs error = %.3e (tol = %.3e).', msg, err, tol));
end