function test_select_primal_dofs()
%TEST_SELECT_PRIMAL_DOFS_CORNERS  Sanity test for corner-based primal selection.
%
% Strengthened checks:
%  (A) Independently recompute expected primal set from geometry and compare sets.
%  (B) Verify glob_c / glob_delta form a disjoint partition of iface.glob.
%  (C) Verify hat_c / hat_delta are consistent with iface.glob2hat mapping.
%
% Note: This test assumes the thesis model setting used throughout:
% - unit square mesh
% - structured subdomains (nSubX x nSubY)
% - Dirichlet at x=0 and x=1 so those nodes are not part of free/interface DOFs
% - primal set = decomposition corner nodes (on both x- and y-partition lines)

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/coarse');

  % Keep mesh divisible by typical nSubX,nSubY so partition lines fall on mesh nodes.
  n = 12;

  % Test several structured decompositions (avoid being correct-by-accident for one case).
  configs = [
    2, 2;
    3, 2;
    4, 3
  ];

  tol = 1e-12;  % must match the fixed tol used inside select_primal_dofs.m

  for ci = 1:size(configs,1)
    nSubX = configs(ci,1);
    nSubY = configs(ci,2);

    [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
    [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY); %#ok<ASGLU>
    [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>

    primal = select_primal_dofs(p, ddm, iface);

    % ---------- (A) Independently compute expected primal set ----------
    xLines = (0:nSubX) / nSubX;
    yLines = (0:nSubY) / nSubY;

    iface_glob = iface.glob(:);
    expected_mask = false(numel(iface_glob), 1);

    for k = 1:numel(iface_glob)
      g = iface_glob(k);
      node = ddm.dof2node(g);
      x = p(node,1);
      y = p(node,2);

      on_x = any(abs(x - xLines) < tol);
      on_y = any(abs(y - yLines) < tol);

      expected_mask(k) = (on_x && on_y);
    end

    expected_glob_c = sort(iface_glob(expected_mask));
    expected_glob_delta = sort(iface_glob(~expected_mask));

    assert(isequal(primal.glob_c(:), expected_glob_c(:)));
    assert(isequal(primal.glob_delta(:), expected_glob_delta(:)));

    % Also keep your original count formula as an extra cheap check:
    expected_count_formula = (nSubX-1) * (nSubY+1);
    assert(primal.nC == expected_count_formula);

    % ---------- (B) Partition invariants ----------
    % Subset property:
    assert(all(ismember(primal.glob_c, iface_glob)));
    assert(all(ismember(primal.glob_delta, iface_glob)));

    % Disjointness:
    assert(isempty(intersect(primal.glob_c, primal.glob_delta)));

    % Union equals the interface set:
    union_sorted = sort([primal.glob_c(:); primal.glob_delta(:)]);
    assert(isequal(union_sorted, sort(iface_glob)));

    % ---------- (C) Hat mapping consistency ----------
    assert(isequal(primal.hat_c(:), iface.glob2hat(primal.glob_c(:))));
    assert(isequal(primal.hat_delta(:), iface.glob2hat(primal.glob_delta(:))));

    fprintf('PASS: test_select_primal_dofs_corners (n=%d, %dx%d subdomains, nC=%d)\n', ...
            n, nSubX, nSubY, primal.nC);
  end
end