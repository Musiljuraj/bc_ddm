function test_build_primal_maps()
%TEST_BUILD_PRIMAL_MAPS  Sanity checks for primal/delta bookkeeping.
%
% This test checks:
%  - glob-space coarse numbering (glob2c) and consistency with select_primal_dofs outputs
%  - hat-space coarse numbering (hat2c) and consistency with select_primal_dofs outputs
%  - per-subdomain local splits in gamma_glob (idx_c/idx_d, glob_c/glob_d, c_ids)
%  - per-subdomain product splits (prod_idx_c/prod_idx_d) and alignment with local ordering
%  - global product-space masks (prod_is_c/prod_is_d) and delta compression maps (delta_idx/prod2delta)

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/interface');
  %addpath('../../src/ddm/coarse');

  % Run a small sweep to avoid overfitting to one decomposition.
  configs = [ ...
    12, 3, 2; ...
    10, 2, 2  ...
  ];

  for c = 1:size(configs,1)
    n     = configs(c,1);
    nSubX = configs(c,2);
    nSubY = configs(c,3);

    [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
    [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
    [sub, iface] = identify_interface_dofs(sub, ddm);
    [sub, prod]  = build_product_interface(sub, iface);

    primal = select_primal_dofs(p, ddm, iface);
    [sub, primal] = build_primal_maps(sub, ddm, iface, prod, primal);

    % ----------------------------
    % Global checks: glob2c + hat2c
    % ----------------------------
    glob_c = primal.glob_c(:);
    nC = numel(glob_c);
    assert(primal.nC == nC);

    % On primal dofs: glob2c must map to 1..nC in the same order as glob_c.
    if nC > 0
      assert(isequal(primal.glob2c(glob_c), (1:nC).'));
    end

    % On explicit delta set from selector: must map to 0.
    if isfield(primal,'glob_delta') && ~isempty(primal.glob_delta)
      assert(all(primal.glob2c(primal.glob_delta(:)) == 0));
    end

    % Hat-side: selector provides hat_c/hat_delta; build_primal_maps must agree.
    if isfield(primal,'hat_c') && ~isempty(primal.hat_c)
      assert(isequal(primal.hat2c(primal.hat_c(:)), (1:nC).'));
    end
    if isfield(primal,'hat_delta') && ~isempty(primal.hat_delta)
      assert(all(primal.hat2c(primal.hat_delta(:)) == 0));
    end

    % Product-side global consistency: prod_is_c should match (hat2c ∘ prod2hat) > 0.
    mapped = primal.hat2c(prod.prod2hat) > 0;
    assert(all(mapped == primal.prod_is_c));
    assert(all(primal.prod_is_d == ~primal.prod_is_c));

    % ----------------------------
    % Per-subdomain checks
    % ----------------------------
    for i = 1:numel(sub)
      % Ensure gamma_glob exists (build_primal_maps normalizes glob_G -> gamma_glob).
      assert(isfield(sub(i),'gamma_glob'));
      gG = sub(i).gamma_glob(:);
      nG = numel(gG);

      % Ensure prod_idx exists and matches length.
      assert(isfield(sub(i),'prod_idx'));
      pidx = sub(i).prod_idx(:);
      assert(numel(pidx) == nG);

      idxc = sub(i).idx_c(:);
      idxd = sub(i).idx_d(:);

      if nG == 0
        % All split fields should exist and be empty.
        assert(isempty(idxc) && isempty(idxd));
        assert(isempty(sub(i).glob_c) && isempty(sub(i).glob_d));
        assert(isempty(sub(i).c_ids));
        assert(isempty(sub(i).prod_idx_c) && isempty(sub(i).prod_idx_d));
        continue;
      end

      % idxc/idxd must partition 1..nG (cover all, disjoint).
      merged_local = sort([idxc; idxd]);
      assert(isequal(merged_local, (1:nG).'));
      assert(isempty(intersect(idxc, idxd)));

      % Local gamma_glob split must agree with stored glob_c/glob_d.
      assert(isequal(gG(idxc), sub(i).glob_c(:)));
      assert(isequal(gG(idxd), sub(i).glob_d(:)));

      % Coarse ids must match global map on the subdomain's primal globals.
      assert(isequal(sub(i).c_ids(:), primal.glob2c(sub(i).glob_c(:))));

      % Product split: prod_idx_c/prod_idx_d must partition prod_idx and align with local ordering.
      idx_prod  = pidx;
      idx_prod_c = sub(i).prod_idx_c(:);
      idx_prod_d = sub(i).prod_idx_d(:);

      merged_prod = sort([idx_prod_c; idx_prod_d]);
      assert(isequal(merged_prod, sort(idx_prod)));
      assert(isempty(intersect(idx_prod_c, idx_prod_d)));

      % Alignment: selecting by local idxc/idxd must reproduce prod_idx_c/prod_idx_d.
      assert(isequal(pidx(idxc), idx_prod_c));
      assert(isequal(pidx(idxd), idx_prod_d));

      % Global product flags must match per-subdomain product splits.
      assert(all(primal.prod_is_c(idx_prod_c)));
      assert(all(~primal.prod_is_c(idx_prod_d)));
    end

    % ----------------------------
    % Delta compression map checks
    % ----------------------------
    assert(isequal(primal.delta_idx(:), find(~primal.prod_is_c)));
    assert(primal.nDeltaProd == numel(primal.delta_idx));

    % prod2delta must be 0 on primal product entries.
    assert(all(primal.prod2delta(primal.prod_is_c) == 0));

    % prod2delta must enumerate delta entries as 1..nDeltaProd in the order delta_idx.
    if primal.nDeltaProd > 0
      assert(isequal(primal.prod2delta(primal.delta_idx), (1:primal.nDeltaProd).'));
    end

    % Count consistency: (#delta) + (#primal) = nProd.
    assert(primal.nDeltaProd + nnz(primal.prod_is_c) == prod.nProd);

    fprintf('PASS: test_build_primal_maps (n=%d, %dx%d subdomains)\n', n, nSubX, nSubY);
  end
end