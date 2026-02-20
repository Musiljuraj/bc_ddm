% test_build_product_interface.m
%
% Unit test for build_product_interface.m (standalone).
% Focus: indexing + mapping correctness for product interface space bookkeeping.
%
% This test is intentionally independent of R/B/primal logic.

function test_build_product_interface()
  fprintf('Running test_build_product_interface...\n');

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/fem/elements');
  %addpath('../../src/fem/assembly');
  %addpath('../../src/fem/bc');
  %addpath('../../src/ddm/local');
  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/interface');

  %------------------------------------------------------------
  % Deterministic RNG (Octave-compatible)
  %------------------------------------------------------------
  try
    rng(1);
  catch
    rand('state', 1);
  end

  %------------------------------------------------------------
  % SECTION A: Synthetic minimal-but-nontrivial interface setup
  % (Validates build_product_interface in isolation)
  %------------------------------------------------------------
  %
  % Global free DOF ids (example): 1..6
  % Interface global DOFs: 2,3,4,5
  %
  % Assembled (hat) indexing:
  %   glob 2 -> hat 1
  %   glob 3 -> hat 2
  %   glob 4 -> hat 3
  %   glob 5 -> hat 4
  %
  % Subdomain interface sets (glob_G) with sharing:
  %   sub1: [2,3]
  %   sub2: [3,4]
  %   sub3: [2,5]
  %   sub4: [4,5]
  %
  % Expected multiplicities in product space:
  %   hat1 (glob2): 2 copies (sub1, sub3)
  %   hat2 (glob3): 2 copies (sub1, sub2)
  %   hat3 (glob4): 2 copies (sub2, sub4)
  %   hat4 (glob5): 2 copies (sub3, sub4)
  %
  maxGlob = 6;
  iface = struct();
  iface.nHat = 4;
  iface.glob2hat = zeros(maxGlob, 1);
  iface.glob2hat(2) = 1;
  iface.glob2hat(3) = 2;
  iface.glob2hat(4) = 3;
  iface.glob2hat(5) = 4;

  sub = repmat(struct('glob_G', []), 4, 1);
  sub(1).glob_G = [2,3];
  sub(2).glob_G = [3,4];
  sub(3).glob_G = [2,5];
  sub(4).glob_G = [4,5];

  [sub_out, prod] = build_product_interface(sub, iface);

  %------------------------------------------------------------
  % (1) DIMENSION CONSISTENCY
  %------------------------------------------------------------
  expected_nHat  = iface.nHat;
  expected_nProd = numel(sub(1).glob_G) + numel(sub(2).glob_G) + numel(sub(3).glob_G) + numel(sub(4).glob_G);

  assert(prod.nHat == expected_nHat,  sprintf('nHat mismatch: got %d expected %d', prod.nHat, expected_nHat));
  assert(prod.nProd == expected_nProd, sprintf('nProd mismatch: got %d expected %d', prod.nProd, expected_nProd));

  % Also verify nProd equals sum_i |glob_G|
  sum_local = 0;
  for i = 1:numel(sub_out)
    sum_local = sum_local + numel(sub_out(i).gamma_glob);
  end
  assert(sum_local == prod.nProd, sprintf('Sum of local gamma sizes != nProd: %d vs %d', sum_local, prod.nProd));

  %------------------------------------------------------------
  % (2) INDEX VALIDITY + CONTIGUITY OF PRODUCT INDICES
  %------------------------------------------------------------
  % Every product index 1..nProd should appear exactly once across sub(i).prod_idx
  all_idx = [];
  for i = 1:numel(sub_out)
    all_idx = [all_idx; sub_out(i).prod_idx(:)];
  end
  all_idx_sorted = sort(all_idx);
  assert(numel(all_idx_sorted) == prod.nProd, 'Total number of collected prod_idx entries != nProd');
  assert(all(all_idx_sorted == (1:prod.nProd).'), 'Product indices are not a permutation of 1..nProd');

  % prod.prod2sub and prod.prod2hat must be fully assigned (no zeros)
  assert(all(prod.prod2sub >= 1 & prod.prod2sub <= numel(sub_out)), 'prod2sub contains out-of-range entries');
  assert(all(prod.prod2hat >= 1 & prod.prod2hat <= prod.nHat),        'prod2hat contains out-of-range entries');

  % Per-subdomain stored fields must align with sizes
  for i = 1:numel(sub_out)
    gi = sub_out(i).gamma_glob;
    hi = sub_out(i).gamma_hat;
    pi = sub_out(i).prod_idx;
    assert(numel(gi) == numel(hi), sprintf('sub(%d): gamma_glob and gamma_hat size mismatch', i));
    assert(numel(gi) == numel(pi), sprintf('sub(%d): gamma_glob and prod_idx size mismatch', i));

    % gamma_hat should equal iface.glob2hat(gamma_glob)
    if ~isempty(gi)
      mapped = iface.glob2hat(gi(:));
      assert(all(mapped == hi(:)), sprintf('sub(%d): gamma_hat != glob2hat(gamma_glob)', i));
    end

    % prod.prod2sub/prod.prod2hat should agree with sub’s local bookkeeping
    if ~isempty(pi)
      assert(all(prod.prod2sub(pi) == i), sprintf('sub(%d): prod2sub mismatch on prod_idx', i));
      assert(all(prod.prod2hat(pi) == hi(:)), sprintf('sub(%d): prod2hat mismatch on prod_idx', i));
    end
  end

  %------------------------------------------------------------
  % (3) MULTIPLICITY CONSISTENCY
  %------------------------------------------------------------
  % Multiplicity by counting prod.prod2hat
  mult_count = zeros(prod.nHat, 1);
  for j = 1:prod.nProd
    mult_count(prod.prod2hat(j)) = mult_count(prod.prod2hat(j)) + 1;
  end

  % Multiplicity by hat2prod list lengths must match
  for hk = 1:prod.nHat
    list = prod.hat2prod{hk};
    assert(numel(list) == mult_count(hk), sprintf('hat2prod multiplicity mismatch for hk=%d', hk));
    if ~isempty(list)
      % Each entry in hat2prod{hk} must map back to hk via prod2hat
      assert(all(prod.prod2hat(list) == hk), sprintf('hat2prod contains wrong indices for hk=%d', hk));
    end
  end

  % Expected multiplicities in this synthetic setup are all 2
  assert(all(mult_count == 2), sprintf('Expected all multiplicities 2, got [%s]', num2str(mult_count.')));

  %------------------------------------------------------------
  % (4) STRONG MAPPING CONSISTENCY CHECK (manual gather)
  %------------------------------------------------------------
  u_hat  = rand(prod.nHat, 1);
  u_prod = zeros(prod.nProd, 1);

  % Manual gather using prod.prod2hat:
  for j = 1:prod.nProd
    u_prod(j) = u_hat(prod.prod2hat(j));
  end

  % Verify consistency via hat2prod lists:
  for hk = 1:prod.nHat
    list = prod.hat2prod{hk};
    for r = 1:numel(list)
      j = list(r);
      assert(u_prod(j) == u_hat(hk), sprintf('Gather mismatch: u_prod(%d) != u_hat(%d)', j, hk));
    end
  end

  fprintf('  Synthetic-case checks: PASSED\n');

  %------------------------------------------------------------
  % SECTION B (optional): If richer iface fields exist in your pipeline,
  % check multiplicities against iface.dof2sub + iface.hat2glob.
  %
  % This block does NOT fail the test if those fields are absent.
  % It adds a stronger check when available.
  %------------------------------------------------------------
  if isfield(iface, 'dof2sub') && isfield(iface, 'hat2glob')
    % Here we assume:
    %   iface.hat2glob(hk) gives the global free DOF id for assembled hk,
    %   iface.dof2sub{g} lists subdomains sharing global free DOF g.
    ok = true;
    for hk = 1:prod.nHat
      g = iface.hat2glob(hk);
      if g <= numel(iface.dof2sub) && ~isempty(iface.dof2sub{g})
        expected = numel(iface.dof2sub{g});
        got = numel(prod.hat2prod{hk});
        if expected ~= got
          ok = false;
          fprintf('  NOTE: multiplicity mismatch vs iface.dof2sub at hk=%d (expected %d got %d)\n', hk, expected, got);
        end
      end
    end
    assert(ok, 'Multiplicity check vs iface.dof2sub failed (fields present but inconsistent)');
    fprintf('  Optional dof2sub-based multiplicity check: PASSED\n');
  else
    fprintf('  Optional dof2sub-based multiplicity check: SKIPPED (iface lacks dof2sub/hat2glob)\n');
  end

  fprintf('test_build_product_interface: ALL PASSED\n');
end

% Allow running as a script:
test_build_product_interface();
