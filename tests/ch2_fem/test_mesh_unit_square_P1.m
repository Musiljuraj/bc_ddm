function test_mesh_unit_square_P1()
%TEST_MESH_UNIT_SQUARE_P1 Basic sanity checks for the mesh generator.

  setup_paths();
  %addpath('../../src/fem/mesh');

  n = 4;
  [p, t, bnd] = mesh_unit_square_P1(n);

  np = n + 1;
  N_expected  = np*np;
  Ne_expected = 2*n*n;

  assert(size(p,1) == N_expected);
  assert(size(p,2) == 2);
  assert(size(t,1) == Ne_expected);
  assert(size(t,2) == 3);

  % Check coordinate bounds
  assert(all(p(:,1) >= -1e-14 & p(:,1) <= 1+1e-14));
  assert(all(p(:,2) >= -1e-14 & p(:,2) <= 1+1e-14));

  % Check Dirichlet nodes are exactly left/right edges: count = 2*np)
  dirichlet = bnd.dirichlet_nodes(:);
  assert(numel(dirichlet) == 2*np);

  % Check Neumann edges count: bottom n + top n = 2n
  assert(size(bnd.neumann_edges,1) == 2*n);
  assert(size(bnd.neumann_edges,2) == 2);

  % Check a couple of corner nodes are in Dirichlet set
  % With left-to-right, top-to-bottom
  % bottom-left is node 1, bottom-right is node np, top-left is node 1+(np-1)*np, top-right is np*np
  bl = 1;
  br = np;
  tl = 1 + (np-1)*np;
  tr = np*np;
  assert(any(dirichlet == bl));
  assert(any(dirichlet == br));
  assert(any(dirichlet == tl));
  assert(any(dirichlet == tr));

  fprintf('test_mesh_unit_square_P1: OK (n=%d)\n', n);
end
