function test_assemble_stiffness_P1()
%TEST_ASSEMBLE_STIFFNESS_P1 Basic sanity checks for global stiffness assembly.

  setup_paths();
  %addpath('../../src/fem/assembly');

  tol = 1e-12;

  n = 4;
  [p, t, bnd] = mesh_unit_square_P1(n);
  K = assemble_stiffness_P1(p, t);

  N = size(p,1);  % number of global degrees of freedom

  % size + sparsity type
  assert(all(size(K) == [N, N]), 'K has wrong size.');
  assert(issparse(K), 'K should be sparse.');

  % symmetry
  assert(norm(K - K.', 'fro') < tol, 'K is not symmetric (unexpected).');

  % constant nullspace property (for the raw stiffness matrix, before Dirichlet nodes removal)
  r = K * ones(N,1);
  assert(norm(r, 2) < tol, 'K*ones != 0 (constant nullspace violated).');

  % diagonal should be nonnegative
  assert(all(diag(K) >= -tol), 'K has significantly negative diagonal entries.');

  fprintf('PASS: test_assemble_stiffness_P1 (n=%d)\n', n);
end
