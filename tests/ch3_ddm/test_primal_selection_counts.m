function test_primal_selection_counts()
%TEST_PRIMAL_SELECTION_COUNTS  Check primal DOF selection for structured decomposition corners (Chapter 3.4.2).

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/coarse');

  n = 12;
  nSubX = 3;
  nSubY = 2;

  [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY); %#ok<ASGLU>
  [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>

  primal = select_primal_dofs(p, ddm, iface);

  % For this thesis model: Dirichlet at x=0 and x=1 => no free nodes there.
  % Decomposition-corner nodes are at x = k/nSubX with k=1..nSubX-1,
  % and y = l/nSubY with l=0..nSubY. Hence expected count:
  expected = (nSubX-1) * (nSubY+1);

  assert(primal.nC == expected);

  fprintf('PASS: test_primal_selection_counts (n=%d, %dx%d subdomains, nC=%d)\n', n, nSubX, nSubY, primal.nC);
end
