function test_build_assembly_operator_R()
%TEST_BUILD_ASSEMBLY_OPERATOR_R  Unit tests for build_assembly_operator_R (assembly/distribution operator).
%
% Scope:
%   - Construction + basic algebraic meaning of R.
%   - Tiny synthetic deterministic prod struct (no mesh pipeline).
%
% Covers:
%   - output shape/sparsity and 0/1 structure
%   - mapping check: R*u_hat == u_hat(prod2hat)
%   - multiplicity checks: sum(R,1) and diag(R'*R)
%   - basic input validation (nargin, missing fields, bad indices)

  ensure_project_paths_();

  % -------------------------
  % A) Synthetic configuration
  % -------------------------
  prod = struct();
  prod.nHat = 5;
  prod.prod2hat = [1; 2; 2; 4; 5; 1; 3];
  prod.nProd = numel(prod.prod2hat);

  R = build_assembly_operator_R(prod);

  % -------------------------
  % Structural properties
  % -------------------------
  assert(issparse(R), 'R must be sparse.');
  assert(isequal(size(R), [prod.nProd, prod.nHat]), 'R has wrong size.');

  % Exactly one nonzero per row (=> nnz == nProd AND each row has exactly one nz)
  assert(nnz(R) == prod.nProd, 'R must have exactly one nonzero per row (nnz == nProd).');

  row_sums = full(sum(R, 2));
  assert(all(row_sums == 1), 'Each row of R must sum to 1 (one "1" per row).');

  row_nnz = full(sum(spones(R), 2));
  assert(all(row_nnz == 1), 'Each row of R must contain exactly one nonzero.');

  nz = nonzeros(R);
  assert(~isempty(nz) && all(nz == 1), 'All nonzeros of R must be exactly 1.');
  assert(all(isfinite(nz)), 'R contains NaN/Inf (unexpected).');

  % -------------------------
  % Strong mapping check
  % -------------------------
  u_hat = (1:prod.nHat).';
  u_prod_expected = u_hat(prod.prod2hat);
  u_prod_actual   = R * u_hat;
  assert_close_(u_prod_actual, u_prod_expected, 0, 0, ...
    'R does not implement u_prod = u_hat(prod2hat).');

  % -------------------------
  % Multiplicity checks
  % -------------------------
  % Column sums equal number of product copies mapping to each hat DOF
  counts = accumarray(prod.prod2hat, 1, [prod.nHat, 1]);
  col_sums = full(sum(R, 1)).';
  assert(isequal(col_sums, counts), ...
    'Column sums of R do not match multiplicities implied by prod2hat.');

  % Equivalent invariant: R''*R is diagonal with these multiplicities
  C = full(R' * R);
  assert(all(all(C - diag(diag(C)) == 0)), 'R''*R must have zero off-diagonals.');
  assert(isequal(diag(C), counts), 'diag(R''*R) must equal multiplicities.');

  % --------------------------------
  % B) Basic input validation checks
  % --------------------------------
  assert_throws_(@() build_assembly_operator_R(),     'expected error for nargin=0');
  assert_throws_(@() build_assembly_operator_R(1,2), 'expected error for nargin=2');
  assert_throws_(@() build_assembly_operator_R(123), 'expected error for non-struct input');

  bad = rmfield(prod, 'prod2hat');
  assert_throws_(@() build_assembly_operator_R(bad), 'expected error for missing prod2hat');

  bad = prod;
  bad.prod2hat = bad.prod2hat(1:end-1);
  assert_throws_(@() build_assembly_operator_R(bad), 'expected error for prod2hat length mismatch');

  bad = prod;
  bad.prod2hat(3) = prod.nHat + 1;
  assert_throws_(@() build_assembly_operator_R(bad), 'expected error for out-of-range prod2hat index');

  fprintf('PASS: test_build_assembly_operator_R\n');
end

% =========================
% Helpers (same style as test_apply_local_schur.m)
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

function assert_close_(a, b, rtol, atol, msg)
  if nargin < 3, rtol = 0; end
  if nargin < 4, atol = 0; end
  if nargin < 5, msg = 'Values not close.'; end

  da = full(a); db = full(b);
  err = norm(da - db, inf);
  denom = max(1, norm(db, inf));
  ok = (err <= atol + rtol * denom);

  assert(ok, sprintf('%s (err=%g, denom=%g, rtol=%g, atol=%g)', msg, err, denom, rtol, atol));
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