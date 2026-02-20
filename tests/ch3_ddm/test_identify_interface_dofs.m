function test_identify_interface_dofs()
%TEST_IDENTIFY_INTERFACE_DOFS  Sanity checks for interior/interface split (Chapter 3.1.3).

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/ddm/partition');

  n = 8;
  nSubX = 2;
  nSubY = 2;

  [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);

  [sub, iface] = identify_interface_dofs(sub, ddm);

  % Global accounting: sum of multiplicities equals sum of local sizes.
  sum_local = 0;
  for i = 1:numel(sub)
    sum_local = sum_local + numel(sub(i).loc2glob);
  end
  assert(abs(sum(iface.counts) - sum_local) < 1e-12);

  % Interface DOFs have multiplicity >= 2.
  assert(all(iface.counts(iface.glob) >= 2));

  % Per subdomain: glob_G are exactly those with counts>1.
  for i = 1:numel(sub)
    g = sub(i).loc2glob(:);
    maskG = iface.counts(g) > 1;

    assert(isequal(sub(i).glob_G(:), g(maskG)));
    assert(isequal(sub(i).glob_I(:), g(~maskG)));
  end

  fprintf('PASS: test_identify_interface_dofs (n=%d, nSubX=%d, nSubY=%d)\n', n, nSubX, nSubY);
end
