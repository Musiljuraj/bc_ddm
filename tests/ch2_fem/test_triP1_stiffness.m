function test_triP1_stiffness()
%TEST_TRIP1_STIFFNESS  Sanity/unit tests for triP1_stiffness (P1 triangle stiffness).
%
% This test suite verifies:
%   1) Correct stiffness matrix, area, basis gradients on the reference triangle, symmetry and constant-nullspace property
%   2) Invariance of Ke under vertex orientation change (CCW vs CW ordering).
%   3) Proper error handling for degenerate triangles (zero/near-zero area).
%
% Usage:
%   test_triP1_stiffness()

  setup_paths();
  %addpath('../../src/fem/elements');
  
  tol = 1e-12;

  % -----------------------------------------------------------------------
  % Test 1: Reference triangle (0,0)-(1,0)-(0,1), CCW ordering
  % -----------------------------------------------------------------------
  xy = [0 0;
        1 0;
        0 1];

  [Ke, area, gradphi] = triP1_stiffness(xy);

  % Expected area: 1/2
  assert(abs(area - 0.5) < tol, 'Area mismatch on reference triangle.');

  % Expected gradients on reference triangle:
  %   phi1 = 1 - x - y  -> grad = [-1; -1]
  %   phi2 = x          -> grad = [ 1;  0]
  %   phi3 = y          -> grad = [ 0;  1]
  grad_expected = [-1  1  0;
                   -1  0  1];
  assert(norm(gradphi - grad_expected, 'fro') < tol, 'gradphi mismatch on reference triangle.');

  % Expected stiffness matrix on reference triangle:
  Ke_expected = [ 1.0  -0.5  -0.5;
                 -0.5   0.5   0.0;
                 -0.5   0.0   0.5];
  assert(norm(Ke - Ke_expected, 'fro') < tol, 'Ke mismatch on reference triangle.');

  % Symmetry check
  assert(norm(Ke - Ke.', 'fro') < tol, 'Ke is not symmetric.');

  % Constant-nullspace check: constants must be in a kernel of the Laplacian
  assert(norm(Ke * ones(3,1), 2) < tol, 'Ke*ones != 0 (constant nullspace violated).');

  % -----------------------------------------------------------------------
  % Test 2: Orientation invariance (same triangle, CW ordering)
  % Ke should remain identical; gradphi may differ by sign due to signed area.
  % -----------------------------------------------------------------------
  xy_cw = [0 0;
           0 1;
           1 0];

  [Ke2, area2, gradphi2] = triP1_stiffness(xy_cw);

  assert(abs(area2 - 0.5) < tol, 'Area mismatch on CW ordering.');
  assert(norm(Ke2 - Ke_expected, 'fro') < tol, 'Ke depends on vertex orientation (should not).');

  % gradphi sign convention depends on ordering; only dot-products are required
  assert(norm((Ke2 - Ke2.'), 'fro') < tol, 'Ke2 is not symmetric.');
  assert(norm(Ke2 * ones(3,1), 2) < tol, 'Ke2*ones != 0 (constant nullspace violated).');

  % -----------------------------------------------------------------------
  % Test 3: Degenerate triangle should raise an error
  % -----------------------------------------------------------------------
  xy_deg = [0 0;
            1 0;
            2 0];  % collinear points => zero area

  did_error = false;
  try
    triP1_stiffness(xy_deg);
  catch
    did_error = true;
  end
  assert(did_error, 'Degenerate triangle did not raise an error.');

  fprintf('PASS: test_triP1_stiffness\n');
end
