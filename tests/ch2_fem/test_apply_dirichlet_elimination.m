function test_apply_dirichlet_elimination()
%TEST_APPLY_DIRICHLET_ELIMINATION  Basic unit tests for apply_dirichlet_elimination.m
%
% THEORY BEHIND THESE TESTS
% ------------------------
% Applying homogeneous Dirichlet BCs by elimination corresponds to restricting
% the assembled linear system K u = F to the free degrees of freedom (DOFs).
% If D is the set of Dirichlet nodes and F is the set of free nodes, then:
%
%   Kff = K(F,F),   Ff = F(F),
%
% and the reduced system is:
%   Kff * u_F = Ff
%
% These tests verify:
%  1) Correct index set of free nodes
%  2) Correct sizes of Kff and Ff
%  3) Correct extracted submatrix/subvector
%  4) Robustness to unsorted and duplicate dirichlet_nodes
%  5) Proper error handling for invalid indices

  setup_paths();
  %addpath('../../src/fem/elements');
  %addpath('../../src/fem/assembly');
  %addpath('../../src/fem/mesh');
  %addpath('../../src/fem/bc');
  
  tol = 1e-12;


  % ------------------------------------------------------------
  % Test 1: Correct extraction on a small hand-constructed system
  % ------------------------------------------------------------
  % Example symmetric matrix (size 5x5) and RHS:
  K = [ 4 -1  0  0  0;
       -1  4 -1  0  0;
        0 -1  4 -1  0;
        0  0 -1  4 -1;
        0  0  0 -1  4 ];

  F = [10; 20; 30; 40; 50];

  % Dirichlet nodes: choose {2,5}
  dirichlet_nodes = [2; 5];

  [Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, F, dirichlet_nodes);

  % Free nodes should be {1,3,4}
  free_expected = [1; 3; 4];
  assert(isequal(free_nodes, free_expected), ...
         'Test 1 failed: free_nodes incorrect.');

  % Sizes: Kff is 3x3, Ff is 3x1
  assert(all(size(Kff) == [3, 3]), ...
         'Test 1 failed: Kff has wrong size.');
  assert(all(size(Ff) == [3, 1]), ...
         'Test 1 failed: Ff has wrong size.');

  % Exact content check:
  Kff_expected = K(free_expected, free_expected);
  Ff_expected  = F(free_expected);

  assert(norm(Kff - Kff_expected, 'fro') < tol, ...
         'Test 1 failed: Kff entries incorrect.');
  assert(norm(Ff - Ff_expected, 2) < tol, ...
         'Test 1 failed: Ff entries incorrect.');

  % ------------------------------------------------------------
  % Test 2: Robustness to duplicates and unsorted dirichlet_nodes
  % ------------------------------------------------------------
  dirichlet_nodes2 = [5; 2; 2; 5];  % same set, different order + duplicates
  [Kff2, Ff2, free_nodes2] = apply_dirichlet_elimination(K, F, dirichlet_nodes2);

  assert(isequal(free_nodes2, free_expected), ...
         'Test 2 failed: free_nodes not robust to duplicates/ordering.');
  assert(norm(Kff2 - Kff_expected, 'fro') < tol, ...
         'Test 2 failed: Kff not robust to duplicates/ordering.');
  assert(norm(Ff2 - Ff_expected, 2) < tol, ...
         'Test 2 failed: Ff not robust to duplicates/ordering.');

  % ------------------------------------------------------------
  % Test 3: Invalid indices must error (index > N)
  % ------------------------------------------------------------
  did_error = false;
  try
    apply_dirichlet_elimination(K, F, [2; 6]); % N=5, so 6 is invalid
  catch
    did_error = true;
  end
  assert(did_error, 'Test 3 failed: invalid index > N did not throw error.');

  % ------------------------------------------------------------
  % Test 4: Invalid indices must error (index < 1)
  % ------------------------------------------------------------
  did_error = false;
  try
    apply_dirichlet_elimination(K, F, [0; 2]);
  catch
    did_error = true;
  end
  assert(did_error, 'Test 4 failed: invalid index < 1 did not throw error.');

  % ------------------------------------------------------------
  % Test 5: Non-integer indices must error
  % ------------------------------------------------------------
  did_error = false;
  try
    apply_dirichlet_elimination(K, F, [2; 3.5]);
  catch
    did_error = true;
  end
  assert(did_error, 'Test 5 failed: non-integer index did not throw error.');

  fprintf('PASS: test_apply_dirichlet_elimination\n');
end
