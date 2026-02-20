function test_assemble_subdomain_matrices_P1()
%TEST_ASSEMBLE_SUBDOMAIN_MATRICES_P1  Verify local subdomain assembly is consistent with global assembly.
%
% Main idea:
%   The global stiffness matrix K and load F are assembled element-wise.
%   If we partition elements into subdomains Ω_i and assemble K^(i), f^(i)
%   from those elements, then summing all subdomain contributions lifted to
%   global FREE DOFs must reproduce the globally assembled reduced system Kff, Ff.

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/fem/elements');
  %addpath('../../src/fem/assembly');
  %addpath('../../src/fem/bc');

  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/local');

  tol = 1e-12;

  n = 6;
  nSubX = 3;
  nSubY = 2;

  [p, t, bnd] = mesh_unit_square_P1(n);

  f = @(x,y) 2.0 + x - 3*y;

  % Global reference
  K = assemble_stiffness_P1(p, t);
  F = assemble_load_P1(p, t, f);
  [Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, full(F), bnd.dirichlet_nodes);

  % Subdomain assembly
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
  [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>
  sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f);

  % Lift + sum into global free DOF matrix/vector
  nFree = ddm.nFree;
  Ksum = sparse(nFree, nFree);
  Fsum = zeros(nFree, 1);

  for i = 1:ddm.Nsub
    g = sub(i).loc2glob(:);
    Ksum(g,g) = Ksum(g,g) + sub(i).K;
    Fsum(g)   = Fsum(g)   + full(sub(i).f);
  end

  % Check ddm free node set matches apply_dirichlet_elimination
  assert(isequal(ddm.dof2node(:), free_nodes(:)), 'Mismatch in free node indexing.');

  % Compare
  assert(norm(Ksum - Kff, 'fro') < tol, 'Sum of local K^(i) does not match global Kff.');
  assert(norm(Fsum - Ff, 2) < tol, 'Sum of local f^(i) does not match global Ff.');

  fprintf('PASS: test_assemble_subdomain_matrices_P1 (n=%d, %dx%d subdomains)\n', n, nSubX, nSubY);
end
