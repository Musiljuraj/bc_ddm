function test_triP1_load()
%TEST_TRIP1_LOAD Unit tests for triP1_load.m (P1 element load vector).
%
% These tests are aligned with the theoretical definition:
%   fe(m) = ∫_T f(x) φ_m(x) dx
% and with the implemented centroid quadrature:
%   fe ≈ f(xc,yc) * |T|/3 * [1;1;1].
    
  setup_paths();
  %addpath('../../src/fem/elements');

  tol = 1e-12;

  % ------------------------------------------------------------
  % Test 1: Reference triangle, constant f (exact for centroid rule)
  % Triangle: (0,0)-(1,0)-(0,1), area = 1/2, centroid = (1/3,1/3)
  % ------------------------------------------------------------
  xy = [0 0; 1 0; 0 1];
  area = 0.5;

  c = 2.5;
  f = @(x,y) c;

  fe = triP1_load(xy, f);

  % Expected: each entry = c * area / 3
  fe_expected = (c * area / 3) * ones(3,1);

  assert(norm(fe - fe_expected, 2) < tol, ...
         'Test 1 failed: constant f on reference triangle.');

  % Sum of entries should equal integral of f over the triangle: c*area
  assert(abs(sum(fe) - c*area) < tol, ...
         'Test 1 failed: sum(fe) does not equal c*area.');

  % ------------------------------------------------------------
  % Test 2: Orientation invariance (CCW vs CW order)
  % ------------------------------------------------------------
  xy_ccw = [0 0; 1 0; 0 1];
  xy_cw  = [0 0; 0 1; 1 0];

  fe_ccw = triP1_load(xy_ccw, f);
  fe_cw  = triP1_load(xy_cw,  f);

  assert(norm(fe_ccw - fe_cw, 2) < tol, ...
         'Test 2 failed: orientation invariance violated.');

  % ------------------------------------------------------------
  % Test 3: Degenerate triangle must throw an error
  % ------------------------------------------------------------
  xy_deg = [0 0; 1 0; 2 0];  % collinear points => area = 0

  did_error = false;
  try
    triP1_load(xy_deg, f);
  catch
    did_error = true;
  end
  assert(did_error, 'Test 4 failed: degenerate triangle did not error.');

  fprintf('PASS: test_triP1_load\n');
end
