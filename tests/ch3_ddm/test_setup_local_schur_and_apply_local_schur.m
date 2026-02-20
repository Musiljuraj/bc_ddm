function test_setup_local_schur_and_apply_local_schur()
%TEST_SETUP_LOCAL_SCHUR_AND_APPLY_LOCAL_SCHUR  Validate Schur complement logic.
%
% Verifies on at least one subdomain:
%   - g^(i) matches the theoretical formula,
%   - matrix-free apply_local_schur matches explicit Schur complement S^(i).

  addpath('../../src/fem/mesh');
  addpath('../../src/fem/elements');

  addpath('../../src/ddm/partition');
  addpath('../../src/ddm/local');

  tol = 1e-10;

  n = 6;
  nSubX = 3;
  nSubY = 2;

  [p, t, bnd] = mesh_unit_square_P1(n);
  f = @(x,y) 1 + x + y;

  [sub, ddm]   = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
  [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>

  sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f);
  sub = extract_subdomain_blocks(sub);

  opts = struct();
  opts.assemble_S = true;    % so we can compare to explicit S
  sub = setup_local_schur(sub, opts);

  % Choose one subdomain that has interior DOFs
  pick = 0;
  for i = 1:ddm.Nsub
    if ~isempty(sub(i).K_II)
      pick = i;
      break;
    end
  end
  assert(pick ~= 0, 'No subdomain with interior DOFs found (unexpected for this test).');

  si = sub(pick);

  % --- test g formula against direct backslash
  g_ref = si.f_g - si.K_gI * (si.K_II \ si.f_I);
  assert(norm(full(si.g) - full(g_ref), 2) < tol, 'Reduced RHS g^(i) mismatch.');

  % --- test matrix-free apply against explicit S
  ng = size(si.K_gg,1);
  x = (1:ng).';          % deterministic test vector (no randomness)
  y1 = apply_local_schur(si, x);
  y2 = si.S * x;

  assert(norm(full(y1) - full(y2), 2) < tol, 'apply_local_schur does not match explicit S^(i).');

  fprintf('PASS: test_setup_local_schur_and_apply_local_schur (subdomain=%d)\n', pick);
end
