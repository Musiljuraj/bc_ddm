function test_build_primal_maps()
%TEST_BUILD_PRIMAL_MAPS  Primal/delta bookkeeping: invariants + basic input contract.
%
% % FIXED : Removed invalid assumption "2x1 partition => empty primal set".
% % ADDED : Deterministic nC==0 coverage by injecting empty primal selection into build_primal_maps.
% % FIXED : rmfield negative-tests now remove fields from the whole struct array (Octave struct-array rule).
% % ADDED : Extra hardened-contract tests: duplicate glob_c, non-interface glob_c, prod2hat/glob2hat mismatch.

  ensure_project_paths_();

  configs = [ ...
    4, 2, 2; ...
    6, 3, 2; ...
    4, 2, 1  ...
  ];

  baseline = struct();
  have_baseline = false;

  for c = 1:size(configs,1)
    n     = configs(c,1);
    nSubX = configs(c,2);
    nSubY = configs(c,3);

    [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
    [sub_in, ddm]   = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
    [sub_in, iface] = identify_interface_dofs(sub_in, ddm);
    [sub_in, prod]  = build_product_interface(sub_in, iface);

    primal_in = select_primal_dofs(p, ddm, iface);

    [sub, primal] = build_primal_maps(sub_in, ddm, iface, prod, primal_in);

    local_check_outputs_(sub, ddm, iface, prod, primal);

    % % ADDED: deterministic nC==0 branch coverage by injection (do once)
    if ~have_baseline
      primal0 = primal_in;
      primal0.glob_c = zeros(0,1);
      primal0.hat_c  = zeros(0,1);

      [sub0, primal_out0] = build_primal_maps(sub_in, ddm, iface, prod, primal0);
      local_check_outputs_(sub0, ddm, iface, prod, primal_out0);

      assert(primal_out0.nC == 0, 'Injected-empty: expected primal_out0.nC == 0.');
      assert(isempty(primal_out0.glob_c), 'Injected-empty: expected primal_out0.glob_c empty.');
      assert(all(primal_out0.glob2c == 0), 'Injected-empty: expected glob2c all zeros.');
      assert(all(primal_out0.hat2c  == 0), 'Injected-empty: expected hat2c all zeros.');
      assert(all(primal_out0.prod_is_c == false), 'Injected-empty: expected prod_is_c all false.');
      assert(isequal(primal_out0.delta_idx(:), (1:prod.nProd).'), 'Injected-empty: expected delta_idx == 1:nProd.');
      assert(isequal(primal_out0.prod2delta(:), (1:prod.nProd).'), 'Injected-empty: expected prod2delta == 1:nProd.');
    end

    if ~have_baseline
      baseline.p         = p;
      baseline.sub_in    = sub_in;
      baseline.ddm       = ddm;
      baseline.iface     = iface;
      baseline.prod      = prod;
      baseline.primal_in = primal_in;
      baseline.sub       = sub;
      baseline.primal    = primal;
      have_baseline = true;
    end

    fprintf('PASS: test_build_primal_maps (n=%d, %dx%d subdomains)\n', n, nSubX, nSubY);
  end

  % ----------------------------
  % Signature checks
  % ----------------------------
  local_assert_throws_any_(@() build_primal_maps(), 'nargin=0');
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod), 'nargin=4');
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, baseline.primal_in, 123), 'nargin=6');
  local_assert_throws_any_(@() local_call_too_many_outputs_(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, baseline.primal_in), 'nargout=3');

  % ----------------------------
  % Input-contract checks (struct fields)
  % ----------------------------
  ddm_bad = baseline.ddm;
  if isfield(ddm_bad, 'nFree'); ddm_bad = rmfield(ddm_bad, 'nFree'); end
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, ddm_bad, baseline.iface, baseline.prod, baseline.primal_in), 'ddm missing nFree');

  iface_bad = baseline.iface;
  if isfield(iface_bad, 'nHat'); iface_bad = rmfield(iface_bad, 'nHat'); end
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, iface_bad, baseline.prod, baseline.primal_in), 'iface missing nHat');

  iface_bad2 = baseline.iface;
  if isfield(iface_bad2, 'glob2hat'); iface_bad2 = rmfield(iface_bad2, 'glob2hat'); end
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, iface_bad2, baseline.prod, baseline.primal_in), 'iface missing glob2hat');

  prod_bad = baseline.prod;
  if isfield(prod_bad, 'nProd'); prod_bad = rmfield(prod_bad, 'nProd'); end
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, prod_bad, baseline.primal_in), 'prod missing nProd');

  prod_bad2 = baseline.prod;
  if isfield(prod_bad2, 'prod2hat'); prod_bad2 = rmfield(prod_bad2, 'prod2hat'); end
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, prod_bad2, baseline.primal_in), 'prod missing prod2hat');

  primal_bad = baseline.primal_in;
  if isfield(primal_bad, 'glob_c'); primal_bad = rmfield(primal_bad, 'glob_c'); end
  local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'primal missing glob_c');

  % sub struct: must contain prod_idx and either gamma_glob or glob_G
  sub_bad = baseline.sub_in;
  % % FIXED: rmfield must be applied to the whole struct array; per-element assignment breaks in Octave.
  if isstruct(sub_bad) && isfield(sub_bad, 'prod_idx')
    sub_bad = rmfield(sub_bad, 'prod_idx');
  end
  local_assert_throws_any_(@() build_primal_maps(sub_bad, baseline.ddm, baseline.iface, baseline.prod, baseline.primal_in), 'sub missing prod_idx');

  sub_bad2 = baseline.sub_in;
  % % FIXED: Remove both gamma_glob and glob_G from all elements to keep fields consistent across the struct array.
  if isstruct(sub_bad2) && isfield(sub_bad2, 'gamma_glob')
    sub_bad2 = rmfield(sub_bad2, 'gamma_glob');
  end
  if isstruct(sub_bad2) && isfield(sub_bad2, 'glob_G')
    sub_bad2 = rmfield(sub_bad2, 'glob_G');
  end
  local_assert_throws_any_(@() build_primal_maps(sub_bad2, baseline.ddm, baseline.iface, baseline.prod, baseline.primal_in), 'sub missing gamma_glob and glob_G');

  % prod_idx length must match gamma_glob/glob_G length
  sub_bad3 = baseline.sub_in;
  ii = local_first_nonempty_gamma_or_globG_(sub_bad3);
  if ii > 0
    nG = local_gamma_length_(sub_bad3(ii));
    if nG > 0
      sub_bad3(ii).prod_idx = sub_bad3(ii).prod_idx(1:max(0, nG-1));
      local_assert_throws_any_(@() build_primal_maps(sub_bad3, baseline.ddm, baseline.iface, baseline.prod, baseline.primal_in), sprintf('prod_idx length mismatch at sub(%d)', ii));
    end
  end

  % ----------------------------
  % Input-contract checks (index validity)
  % ----------------------------
  if baseline.ddm.nFree >= 1
    primal_bad = baseline.primal_in;
    primal_bad.glob_c = baseline.ddm.nFree + 1;
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c out of range');

    primal_bad = baseline.primal_in;
    primal_bad.glob_c = 0;
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c contains 0');

    primal_bad = baseline.primal_in;
    primal_bad.glob_c = -1;
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c negative');

    primal_bad = baseline.primal_in;
    primal_bad.glob_c = 1.25;
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c non-integer');

    primal_bad = baseline.primal_in;
    primal_bad.glob_c = NaN;
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c NaN');

    primal_bad = baseline.primal_in;
    primal_bad.glob_c = complex(1, 1);
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c complex');

    primal_bad = baseline.primal_in;
    primal_bad.glob_c = true(3,1);
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c logical wrong length');

    % duplicates in glob_c must error (hardened contract)
    primal_bad = baseline.primal_in;
    primal_bad.glob_c = [1; 1];
    local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c duplicates');

    % glob_c containing non-interface dof must error (glob2hat == 0)
    non_iface = find(baseline.iface.glob2hat(:) == 0, 1);
    if ~isempty(non_iface)
      primal_bad = baseline.primal_in;
      primal_bad.glob_c = non_iface;
      local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'glob_c non-interface');
    end

    % hat_c length mismatch
    primal_bad = baseline.primal_in;
    nC = numel(primal_bad.glob_c(:));
    if nC >= 1
      hat_first = baseline.iface.glob2hat(primal_bad.glob_c(1));
      primal_bad.hat_c = [baseline.iface.glob2hat(primal_bad.glob_c(:)); hat_first];
      local_assert_throws_any_(@() build_primal_maps(baseline.sub_in, baseline.ddm, baseline.iface, baseline.prod, primal_bad), 'hat_c length mismatch');
    end
  end

  % ----------------------------
  % Hardened mapping consistency checks (subdomain tampering)
  % ----------------------------
  sub_bad4 = baseline.sub_in;
  jj = local_first_nonempty_gamma_or_globG_(sub_bad4);
  if jj > 0
    if ~isfield(sub_bad4(jj),'gamma_glob') && isfield(sub_bad4(jj),'glob_G')
      sub_bad4(jj).gamma_glob = sub_bad4(jj).glob_G(:);
    end
    if numel(sub_bad4(jj).prod_idx) >= 2
      tmp = sub_bad4(jj).prod_idx(1);
      sub_bad4(jj).prod_idx(1) = sub_bad4(jj).prod_idx(2);
      sub_bad4(jj).prod_idx(2) = tmp;
      local_assert_throws_any_(@() build_primal_maps(sub_bad4, baseline.ddm, baseline.iface, baseline.prod, baseline.primal_in), sprintf('prod_idx tamper mismatch at sub(%d)', jj));
    end
  end

  sub_bad5 = baseline.sub_in;
  non_iface = find(baseline.iface.glob2hat(:) == 0, 1);
  kk = local_first_nonempty_gamma_or_globG_(sub_bad5);
  if ~isempty(non_iface) && kk > 0
    if ~isfield(sub_bad5(kk),'gamma_glob') && isfield(sub_bad5(kk),'glob_G')
      sub_bad5(kk).gamma_glob = sub_bad5(kk).glob_G(:);
    end
    if numel(sub_bad5(kk).gamma_glob) >= 1
      sub_bad5(kk).gamma_glob(1) = non_iface;
      local_assert_throws_any_(@() build_primal_maps(sub_bad5, baseline.ddm, baseline.iface, baseline.prod, baseline.primal_in), sprintf('gamma_glob non-interface at sub(%d)', kk));
    end
  end
end

% ============================================================
% Helpers
% ============================================================

function local_check_outputs_(sub, ddm, iface, prod, primal)
  local_assert_has_fields_(ddm,  {'nFree'});
  local_assert_has_fields_(iface, {'nHat','glob2hat'});
  local_assert_has_fields_(prod, {'nProd','prod2hat'});

  nFree = ddm.nFree;
  nHat  = iface.nHat;
  nProd = prod.nProd;

  local_assert_has_fields_(primal, {'glob_c','glob2c','hat2c','prod_is_c','prod_is_d','delta_idx','prod2delta','nDeltaProd','nC','hat_c'});

  assert(isvector(primal.glob2c) && numel(primal.glob2c) == nFree, 'primal.glob2c must be length nFree.');
  assert(size(primal.glob2c,2) == 1, 'primal.glob2c must be a column vector.');
  assert(all(isfinite(primal.glob2c(:))), 'primal.glob2c must be finite.');
  assert(all(primal.glob2c(:) == round(primal.glob2c(:))), 'primal.glob2c must be integer-valued.');
  assert(all(primal.glob2c(:) >= 0), 'primal.glob2c must be >= 0.');

  assert(isvector(primal.hat2c) && numel(primal.hat2c) == nHat, 'primal.hat2c must be length nHat.');
  assert(size(primal.hat2c,2) == 1, 'primal.hat2c must be a column vector.');
  assert(all(isfinite(primal.hat2c(:))), 'primal.hat2c must be finite.');
  assert(all(primal.hat2c(:) == round(primal.hat2c(:))), 'primal.hat2c must be integer-valued.');
  assert(all(primal.hat2c(:) >= 0), 'primal.hat2c must be >= 0.');

  glob_c = primal.glob_c(:);
  nC = numel(glob_c);
  assert(primal.nC == nC, 'primal.nC must match numel(primal.glob_c).');

  if nC > 0
    assert(isequal(primal.glob2c(glob_c), (1:nC).'), 'glob2c must map glob_c to 1..nC in-order.');
    assert(all(primal.hat2c(:) <= nC), 'hat2c must be <= nC.');
    assert(isequal(primal.hat_c(:), iface.glob2hat(glob_c)), 'hat_c must equal iface.glob2hat(glob_c).');
  else
    assert(all(primal.glob2c(:) == 0), 'If nC==0, glob2c must be all zeros.');
    assert(all(primal.hat2c(:) == 0), 'If nC==0, hat2c must be all zeros.');
    assert(isempty(primal.hat_c), 'If nC==0, hat_c must be empty.');
  end

  assert(isvector(primal.prod_is_c) && numel(primal.prod_is_c) == nProd, 'prod_is_c must be length nProd.');
  assert(isvector(primal.prod_is_d) && numel(primal.prod_is_d) == nProd, 'prod_is_d must be length nProd.');
  assert(all(primal.prod_is_d(:) == ~primal.prod_is_c(:)), 'prod_is_d must be logical complement of prod_is_c.');

  mapped = primal.hat2c(prod.prod2hat(:)) > 0;
  assert(all(mapped(:) == primal.prod_is_c(:)), 'prod_is_c must match (hat2c(prod2hat)>0).');

  assert(isequal(primal.delta_idx(:), find(~primal.prod_is_c(:))), 'delta_idx must be find(~prod_is_c).');
  assert(primal.nDeltaProd == numel(primal.delta_idx), 'nDeltaProd must match numel(delta_idx).');

  assert(isvector(primal.prod2delta) && numel(primal.prod2delta) == nProd, 'prod2delta must be length nProd.');
  assert(size(primal.prod2delta,2) == 1, 'prod2delta must be a column vector.');
  assert(all(isfinite(primal.prod2delta(:))), 'prod2delta must be finite.');
  assert(all(primal.prod2delta(:) == round(primal.prod2delta(:))), 'prod2delta must be integer-valued.');
  assert(all(primal.prod2delta(primal.prod_is_c) == 0), 'prod2delta must be 0 on primal product entries.');

  if primal.nDeltaProd > 0
    assert(isequal(primal.prod2delta(primal.delta_idx), (1:primal.nDeltaProd).'), 'prod2delta must enumerate delta entries 1..nDeltaProd.');
  else
    assert(all(primal.prod2delta(:) == 0), 'If nDeltaProd==0, prod2delta must be all zeros.');
  end
  assert(primal.nDeltaProd + nnz(primal.prod_is_c) == nProd, 'Count mismatch: nDeltaProd + nnz(prod_is_c) must equal nProd.');

  for i = 1:numel(sub)
    assert(isfield(sub(i),'gamma_glob'), sprintf('sub(%d) must contain gamma_glob after build_primal_maps.', i));
    assert(isfield(sub(i),'prod_idx'), sprintf('sub(%d) must contain prod_idx.', i));

    gG = sub(i).gamma_glob(:);
    nG = numel(gG);
    pidx = sub(i).prod_idx(:);

    assert(numel(pidx) == nG, sprintf('sub(%d) prod_idx length must match gamma_glob.', i));

    local_assert_int_vector_(sub(i).idx_c(:), sprintf('sub(%d).idx_c', i));
    local_assert_int_vector_(sub(i).idx_d(:), sprintf('sub(%d).idx_d', i));

    idxc = sub(i).idx_c(:);
    idxd = sub(i).idx_d(:);

    if nG == 0
      assert(isempty(idxc) && isempty(idxd), sprintf('sub(%d): expected empty idx_c/idx_d for empty gamma_glob.', i));
      assert(isempty(sub(i).glob_c) && isempty(sub(i).glob_d), sprintf('sub(%d): expected empty glob_c/glob_d for empty gamma_glob.', i));
      assert(isempty(sub(i).c_ids), sprintf('sub(%d): expected empty c_ids for empty gamma_glob.', i));
      continue;
    end

    assert(all(idxc >= 1 & idxc <= nG), sprintf('sub(%d).idx_c must be within 1..nG.', i));
    assert(all(idxd >= 1 & idxd <= nG), sprintf('sub(%d).idx_d must be within 1..nG.', i));

    assert(all(diff(idxc) >= 1) || isempty(idxc), sprintf('sub(%d).idx_c must be strictly increasing.', i));
    assert(all(diff(idxd) >= 1) || isempty(idxd), sprintf('sub(%d).idx_d must be strictly increasing.', i));

    merged_local = sort([idxc; idxd]);
    assert(isequal(merged_local, (1:nG).'), sprintf('sub(%d): idx_c/idx_d must partition 1..nG.', i));
    assert(isempty(intersect(idxc, idxd)), sprintf('sub(%d): idx_c and idx_d must be disjoint.', i));

    assert(isequal(gG(idxc), sub(i).glob_c(:)), sprintf('sub(%d): glob_c must equal gamma_glob(idx_c).', i));
    assert(isequal(gG(idxd), sub(i).glob_d(:)), sprintf('sub(%d): glob_d must equal gamma_glob(idx_d).', i));

    assert(isequal(sub(i).c_ids(:), primal.glob2c(sub(i).glob_c(:))), sprintf('sub(%d): c_ids must match glob2c(glob_c).', i));

    idx_prod_c = sub(i).prod_idx_c(:);
    idx_prod_d = sub(i).prod_idx_d(:);

    merged_prod = sort([idx_prod_c; idx_prod_d]);
    assert(isequal(merged_prod, sort(pidx)), sprintf('sub(%d): prod_idx_c/prod_idx_d must partition prod_idx.', i));
    assert(isempty(intersect(idx_prod_c, idx_prod_d)), sprintf('sub(%d): prod_idx_c and prod_idx_d must be disjoint.', i));

    assert(isequal(pidx(idxc), idx_prod_c), sprintf('sub(%d): prod_idx_c must match prod_idx(idx_c).', i));
    assert(isequal(pidx(idxd), idx_prod_d), sprintf('sub(%d): prod_idx_d must match prod_idx(idx_d).', i));

    assert(all(primal.prod_is_c(idx_prod_c)), sprintf('sub(%d): prod_idx_c must be marked primal globally.', i));
    assert(all(~primal.prod_is_c(idx_prod_d)), sprintf('sub(%d): prod_idx_d must be marked delta globally.', i));
  end
end

function local_call_too_many_outputs_(sub, ddm, iface, prod, primal)
  [~, ~, ~] = build_primal_maps(sub, ddm, iface, prod, primal);
end

function ii = local_first_nonempty_gamma_or_globG_(sub)
  ii = 0;
  for k = 1:numel(sub)
    nG = local_gamma_length_(sub(k));
    if nG > 0
      ii = k;
      return;
    end
  end
end

function nG = local_gamma_length_(s)
  nG = 0;
  if isfield(s, 'gamma_glob')
    nG = numel(s.gamma_glob(:));
    return;
  end
  if isfield(s, 'glob_G')
    nG = numel(s.glob_G(:));
    return;
  end
end

function local_assert_has_fields_(s, fields)
  for k = 1:numel(fields)
    f = fields{k};
    assert(isfield(s, f), sprintf('Missing required field: %s', f));
  end
end

function local_assert_int_vector_(v, name)
  assert(isvector(v), sprintf('%s must be a vector.', name));
  if isempty(v)
    return;
  end
  assert(isnumeric(v) && isreal(v), sprintf('%s must be real numeric.', name));
  assert(all(isfinite(v(:))), sprintf('%s must be finite.', name));
  assert(all(v(:) == round(v(:))), sprintf('%s must be integer-valued.', name));
end

function local_assert_throws_any_(fh, label)
  did_throw = false;
  try
    fh();
  catch
    did_throw = true;
  end
  if ~did_throw
    fprintf(2, 'NO ERROR THROWN for bad case: %s\n', label);
  end
  assert(did_throw, sprintf('Expected an error but none was thrown (%s).', label));
end

function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did
    return;
  end

  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    if ~isempty(sp)
      maindir = fileparts(sp);
      rootdir = fileparts(maindir);
      addpath(genpath(fullfile(rootdir, 'tests')));
    end
    did = true;
    return;
  end

  thisdir  = fileparts(mfilename('fullpath'));
  testsdir = fileparts(thisdir);
  rootdir  = fileparts(testsdir);
  maindir  = fullfile(rootdir, 'main');

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, 'tests')));

  did = true;
end