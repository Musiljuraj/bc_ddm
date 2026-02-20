function test_B_times_R_is_zero()
%TEST_B_TIMES_R_IS_ZERO  Verify that range(R) is contained in ker(B), i.e. B*R = 0 (Chapter 3.3.3).

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/interface');

  rand('seed', 1); %#ok<RAND>

  n = 8;
  nSubX = 2;
  nSubY = 2;

  [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY); %#ok<ASGLU>
  [sub, iface] = identify_interface_dofs(sub, ddm);

  [sub, prod] = build_product_interface(sub, iface); %#ok<ASGLU>
  R = build_assembly_operator_R(prod);
  B = build_jump_operator_B(prod);

  nHat = prod.nHat;
  w_hat = rand(nHat, 1);
  w = R * w_hat;

  r = B * w;
  tol = 1e-12;
  assert(norm(full(r), 2) < tol);

  fprintf('PASS: test_B_times_R_is_zero (n=%d, %dx%d subdomains)\n', n, nSubX, nSubY);
end
