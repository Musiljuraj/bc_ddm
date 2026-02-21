function test_triP1_stiffness()
%TEST_TRIP1_STIFFNESS  Stronger unit tests for triP1_stiffness (P1 triangle stiffness).
%
% This suite focuses on:
%   - Signature and basic contract checks
%   - Exact reference-triangle values
%   - Generic (scalene) triangle correctness via an independent reference construction
%   - Permutation consistency (vertex relabeling => row/col permutation)
%   - Rigid-motion + uniform-scaling invariances
%   - Degeneracy handling
%   - % ADDED: stricter input validation expectations (numeric/finite/real/non-logical)
%
% Usage:
%   test_triP1_stiffness()

  setup_paths();

  % -----------------------------------------------------------------------
  % % ADDED: Signature checks (guards accidental signature drift)
  % -----------------------------------------------------------------------
  assert(nargin('triP1_stiffness') == 1, 'triP1_stiffness: expected nargin == 1');
  assert(nargout('triP1_stiffness') == 3, 'triP1_stiffness: expected nargout == 3');

  % -----------------------------------------------------------------------
  % % ADDED: Wrong-argument-count must error (Octave may error before entering body)
  % -----------------------------------------------------------------------
  assert_throws(@() triP1_stiffness());
  assert_throws(@() triP1_stiffness(zeros(3,2), 1));

  % Tolerances (scale-aware)
  rtol = 5e-13;
  atol = 5e-14;

  % -----------------------------------------------------------------------
  % Test 1 (kept, but uses common helpers): Reference triangle, CCW ordering
  % -----------------------------------------------------------------------
  xy_ref = [0 0;
            1 0;
            0 1];

  [Ke, area, gradphi] = triP1_stiffness(xy_ref);

  % Expected area: 1/2
  assert_close(area, 0.5, rtol, atol, 'reference triangle area');

  % Expected gradients on reference triangle:
  %   phi1 = 1 - x - y  -> grad = [-1; -1]
  %   phi2 = x          -> grad = [ 1;  0]
  %   phi3 = y          -> grad = [ 0;  1]
  grad_expected = [-1  1  0;
                   -1  0  1];
  assert_close(gradphi, grad_expected, rtol, atol, 'reference triangle gradphi');

  % Expected stiffness matrix on reference triangle:
  Ke_expected = [ 1.0  -0.5  -0.5;
                 -0.5   0.5   0.0;
                 -0.5   0.0   0.5];
  assert_close(Ke, Ke_expected, rtol, atol, 'reference triangle Ke');

  % Symmetry and constant-nullspace checks
  assert_close(Ke, Ke.', rtol, atol, 'Ke symmetry');
  assert_close(Ke * ones(3,1), zeros(3,1), rtol, atol, 'Ke*ones == 0');

  % -----------------------------------------------------------------------
  % % CHANGED: Orientation/permutation test now uses a *scalene* triangle so
  % it cannot pass accidentally by symmetry.
  % -----------------------------------------------------------------------
  xy = [ 0.10  0.20;
         1.30 -0.10;
         0.40  0.90 ];  % scalene, non-right, non-isosceles

  [Ke0, area0, grad0] = triP1_stiffness(xy);

  % Independent reference: build affine P1 basis by solving interpolation system
  [area_ref, grad_ref, Ke_ref] = reference_quantities(xy);

  assert_close(area0, area_ref, rtol, atol, 'scalene area');
  assert_close(grad0, grad_ref, 50*rtol, 50*atol, 'scalene gradphi (reference solve)'); % slightly looser (solve)
  assert_close(Ke0, Ke_ref, 50*rtol, 50*atol, 'scalene Ke (reference solve)');

  % Additional structural invariants
  assert_close(sum(grad0, 2), [0; 0], rtol, atol, 'sum_i grad(phi_i) == 0');
  assert_close(Ke0 * ones(3,1), zeros(3,1), rtol, atol, 'Ke0*ones == 0');
  assert_close(Ke0, Ke0.', rtol, atol, 'Ke0 symmetry');

  % Positive semidefinite with 1D nullspace (constants)
  v = [1; -2; 1]; % sums to 0, should have positive energy for non-degenerate triangle
  assert((v.' * Ke0 * v) > 0, 'Expected v''*Ke*v > 0 for v with sum(v)=0');

  % -----------------------------------------------------------------------
  % % ADDED: Permutation consistency
  % Reordering vertices should permute basis functions, hence permute rows/cols.
  % -----------------------------------------------------------------------
  perms = {
    [1 2 3]   % identity
    [2 3 1]   % even permutation
    [1 3 2]   % odd permutation (orientation flip)
    [3 2 1]   % odd permutation
  };

  for k = 1:numel(perms)
    p = perms{k};
    [Ke_p, area_p, grad_p] = triP1_stiffness(xy(p,:));

    assert_close(area_p, area0, rtol, atol, sprintf('perm %d: area', k));
    assert_close(Ke_p, Ke0(p,p), 50*rtol, 50*atol, sprintf('perm %d: Ke permutation', k));
    assert_close(grad_p, grad0(:,p), 50*rtol, 50*atol, sprintf('perm %d: gradphi permutation', k));
  end

  % -----------------------------------------------------------------------
  % % ADDED: Translation invariance
  % -----------------------------------------------------------------------
  shift = [3.1, -0.7];
  xy_t = xy + ones(3,1) * shift;
  [Ke_t, area_t, grad_t] = triP1_stiffness(xy_t);

  assert_close(area_t, area0, rtol, atol, 'translation: area');
  assert_close(Ke_t, Ke0, 50*rtol, 50*atol, 'translation: Ke');
  assert_close(grad_t, grad0, 50*rtol, 50*atol, 'translation: gradphi');

  % -----------------------------------------------------------------------
  % % ADDED: Rotation invariance (rigid motion)
  % Ke should be invariant; gradients should rotate with the geometry.
  % -----------------------------------------------------------------------
  theta = 0.37 * pi;
  R = [cos(theta) -sin(theta);
       sin(theta)  cos(theta)];
  xy_r = (R * xy.').';
  [Ke_r, area_r, grad_r] = triP1_stiffness(xy_r);

  assert_close(area_r, area0, 50*rtol, 50*atol, 'rotation: area');
  assert_close(Ke_r, Ke0, 200*rtol, 200*atol, 'rotation: Ke');
  assert_close(grad_r, R * grad0, 200*rtol, 200*atol, 'rotation: gradphi transforms');

  % -----------------------------------------------------------------------
  % % ADDED: Uniform scaling invariance
  % For Laplacian stiffness with P1 elements, uniform scaling leaves Ke invariant.
  % -----------------------------------------------------------------------
  s = 2.7;
  [Ke_s, area_s, grad_s] = triP1_stiffness(s * xy);

  assert_close(area_s, s*s*area0, 50*rtol, 50*atol, 'scaling: area');
  assert_close(Ke_s, Ke0, 200*rtol, 200*atol, 'scaling: Ke invariant');
  assert_close(grad_s, grad0 / s, 200*rtol, 200*atol, 'scaling: gradphi scales as 1/s');

  % -----------------------------------------------------------------------
  % Test 3 (kept): Degenerate triangles must error
  % -----------------------------------------------------------------------
  xy_deg = [0 0;
            1 0;
            2 0];  % collinear => zero area
  assert_throws(@() triP1_stiffness(xy_deg));

  % % ADDED: "near-degenerate" example (triggers current absolute threshold)
  xy_near = [0 0;
             1 0;
             1e-15 1e-15];
  assert_throws(@() triP1_stiffness(xy_near));

  % -----------------------------------------------------------------------
  % % ADDED: Stricter input validation expectations
  % These are intended to catch the same class of issues found earlier in
  % mesh_unit_square_P1.m (char/logical/Inf/NaN silently accepted).
  % If these fail, it indicates triP1_stiffness should be hardened.
  % -----------------------------------------------------------------------
  bad_inputs = {
    zeros(2,2)                      % wrong shape
    zeros(3,3)                      % wrong shape
    zeros(3,2,1)                    % not a 2D matrix
    struct('a', 1)                  % wrong type
    ['1' '2'; '3' '4'; '5' '6']     % char 3x2 (BUG if accepted)
    true(3,2)                       % logical 3x2 (BUG if accepted)
    [0 0; 1 0; NaN 1]               % non-finite
    [0 0; 1 0; Inf 1]               % non-finite
    [0 0; 1 0; 1 1i]                % complex
  };

  for k = 1:numel(bad_inputs)
    xy_bad = bad_inputs{k};
    assert_throws(@() triP1_stiffness(xy_bad));
  end

  fprintf('PASS: test_triP1_stiffness\n');
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function assert_close(A, B, rtol, atol, msg)
  if nargin < 5
    msg = 'values not close';
  end
  d = norm(A(:) - B(:), 2);
  nB = norm(B(:), 2);
  tol = atol + rtol * nB;
  if ~(d <= tol)
    error('assert_close:fail', sprintf('%s (||A-B||=%g, tol=%g)', msg, d, tol));
  end
end

function assert_throws(fhandle)
  did_error = false;
  try
    fhandle();
  catch
    did_error = true;
  end
  assert(did_error, 'Expected an error, but no error was thrown.');
end

function [area, gradphi, Ke] = reference_quantities(xy)
  % Independent construction of P1 basis gradients via interpolation solve.
  % For each i:
  %   phi_i(x,y) = a_i + b_i*x + c_i*y, with phi_i(x_j,y_j)=delta_ij.
  x = xy(:,1);
  y = xy(:,2);

  M = [ones(3,1), x, y];   % 3x3
  C = M \ eye(3);          % columns: [a_i; b_i; c_i]
  gradphi = C(2:3, :);     % 2x3

  doubleA = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
  area = 0.5 * abs(doubleA);

  Ke = area * (gradphi.' * gradphi);
end