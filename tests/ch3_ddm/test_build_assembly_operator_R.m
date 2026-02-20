function test_build_assembly_operator_R()

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/fem/elements');
  %addpath('../../src/fem/assembly');
  %addpath('../../src/fem/bc');
  %addpath('../../src/ddm/local');
  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/interface');


%TEST_BUILD_ASSEMBLY_OPERATOR_R  Unit test for build_assembly_operator_R.m
%
% Scope:
%   - Tests ONLY construction and basic algebraic meaning of R.
%   - Uses a synthetic deterministic prod struct (no mesh pipeline).

  % -------------------------
  % A) Synthetic configuration
  % -------------------------
  prod = struct();
  prod.nHat = 5;
  prod.prod2hat = [1; 2; 2; 4; 5; 1; 3];
  prod.nProd = numel(prod.prod2hat);

  R = build_assembly_operator_R(prod);



  % Structural properties
  assert(isequal(size(R), [prod.nProd, prod.nHat]), 'R has wrong size.');
  assert(nnz(R) == prod.nProd, 'R must have exactly one nonzero per row (nnz == nProd).');

  row_sums = full(sum(R, 2));
  assert(all(row_sums == 1), 'Each row of R must sum to 1 (one "1" per row).');

  row_nnz = full(sum(spones(R), 2));
  assert(all(row_nnz == 1), 'Each row of R must contain exactly one nonzero.');

  assert(all(nonzeros(R) >= 0), 'R must be nonnegative.');

  % Strong mapping check
  u_hat = (1:prod.nHat).';
  u_prod_expected = u_hat(prod.prod2hat);
  u_prod_actual   = R * u_hat;
  assert(norm(u_prod_actual - u_prod_expected, inf) == 0, ...
         'R does not implement u_prod = u_hat(prod2hat).');

  % Column multiplicity check:
  % number of product copies mapping to each hat DOF equals column sum of R
  counts = accumarray(prod.prod2hat, 1, [prod.nHat, 1]);
  col_sums = full(sum(R, 1)).';
  assert(isequal(col_sums, counts), ...
         'Column sums of R do not match multiplicities implied by prod2hat.');

  % --------------------------------
  % B) Basic input validation checks
  % --------------------------------
  assert_throws(@() build_assembly_operator_R(), ...
                'expected input', ...
                'Expected an error when called with no inputs.');

  bad = prod;
  bad.prod2hat = bad.prod2hat(1:end-1);
  assert_throws(@() build_assembly_operator_R(bad), ...
                'wrong length', ...
                'Expected an error for prod2hat length mismatch.');

  bad = prod;
  bad.prod2hat(3) = prod.nHat + 1;
  assert_throws(@() build_assembly_operator_R(bad), ...
                'out-of-range', ...
                'Expected an error for out-of-range prod2hat indices.');

  % Optional final line (consistent with lightweight test scripts)
  disp('test_build_assembly_operator_R: passed');
end

function assert_throws(fhandle, msg_substring, fail_msg)
%ASSERT_THROWS  Assert that fhandle throws and message contains substring.
  did_throw = false;
  try
    fhandle();
  catch err
    did_throw = true;
    % Octave: err.message
    assert(~isempty(strfind(err.message, msg_substring)), ...
           sprintf('Error thrown, but message did not contain "%s". Actual: %s', ...
                   msg_substring, err.message));
  end
  assert(did_throw, fail_msg);
end
