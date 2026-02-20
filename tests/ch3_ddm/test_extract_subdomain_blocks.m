function test_extract_subdomain_blocks()
%TEST_EXTRACT_SUBDOMAIN_BLOCKS  Verify that extracted blocks reconstruct the local matrix.
%
% Since extract_subdomain_blocks is "just slicing", the key failure mode is
% indexing mistakes (wrong orientation / wrong set). This test reconstructs
% K^(i) from blocks and checks equality.

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/fem/elements');

  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/local');

  tol = 1e-12;

  n = 4;
  nSubX = 2;
  nSubY = 2;

  [p, t, bnd] = mesh_unit_square_P1(n);
  [sub, ddm]  = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
  [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>

  f = @(x,y) 1;
  sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f);
  sub = extract_subdomain_blocks(sub);

  for i = 1:ddm.Nsub
    I = sub(i).dofs_I(:);
    G = sub(i).dofs_G(:);

    perm = [I; G];
    Kp = sub(i).K(perm, perm);

    Krec = [sub(i).K_II, sub(i).K_Ig;
            sub(i).K_gI, sub(i).K_gg];

    assert(norm(Kp - Krec, 'fro') < tol, 'Block reconstruction failed on subdomain %d.', i);
  end

  fprintf('PASS: test_extract_subdomain_blocks (n=%d, %dx%d subdomains)\n', n, nSubX, nSubY);
end
