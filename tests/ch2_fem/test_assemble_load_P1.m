function test_assemble_load_P1()
%TEST_ASSEMBLE_LOAD_P1  Basic unit tests for assemble_load_P1.m
%
% THEORY BEHIND THE TESTS
% ----------------------
% Global load vector definition (Galerkin/Ritz):
%   F_i = ∫_Ω φ_i(x) f(x) dx.
%
% Element-wise assembly:
%   F_i = Σ_e ∫_{T_e} φ_i(x) f(x) dx.
%
% With the implemented element quadrature (centroid rule) in triP1_load:
%   fe ≈ f(x_c) * |T|/3 * [1;1;1].
%
% For constant f(x,y) = c, this yields:
%   - sum(fe) = c*|T| for each element (exact under this quadrature)
% Therefore, summing over all elements:
%   sum(F) = Σ_e sum(fe) = c * Σ_e |T_e| = c * |Ω| = c
% because the mesh exactly partitions Ω = [0,1]×[0,1] with total area 1.
%
% This is used as a clean, global consistency check that is independent of
% node ordering and local-to-global indexing details.
  setup_paths();
  %addpath('../../src/fem/elements');
  %addpath('../../src/fem/assembly');
  %addpath('../../src/fem/mesh');
  
  tol = 1e-12;

  % ------------------------------------------------------------
  % Test 1: Dimension check (n=4 => N=(n+1)^2)
  % ------------------------------------------------------------
  n = 4;
  [p, t, bnd] = mesh_unit_square_P1(n); 
  N = size(p,1);

  c = 2.5;
  f = @(x,y) c;

  F = assemble_load_P1(p, t, f);

  assert(all(size(F) == [N, 1]), ...
         'Test 1 failed: F has wrong size (must be N×1).');

  % ------------------------------------------------------------
  % Test 2: Global integral consistency for constant source f=c
  % Expected: sum(F) = ∫_Ω f dx = c * |Ω| = c (since |Ω|=1)
  % ------------------------------------------------------------
  assert(abs(sum(full(F)) - c) < tol, ...
         'Test 2 failed: sum(F) is not equal to integral of f over Ω.');

  % ------------------------------------------------------------
  % Test 3: Zero source term => F must be zero vector
  % ------------------------------------------------------------
  f0 = @(x,y) 0;
  F0 = assemble_load_P1(p, t, f0);

  assert(norm(full(F0), 2) < tol, ...
         'Test 4 failed: zero source does not produce zero load vector.');

  fprintf('PASS: test_assemble_load_P1 (n=%d)\n', n);
end
