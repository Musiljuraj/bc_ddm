function test_build_product_interface()
%TEST_BUILD_PRODUCT_INTERFACE  Unit tests for build_product_interface (product-space assembly).
%
% % CHANGED: deterministic tiny synthetic fixtures (no FEM/partition deps).
% % CHANGED: path-robust helper modeled after test_setup_local_schur.m.
% % ADDED: cover multi-subdomain array incl. one empty glob_G while nHat>0.
% % ADDED: cover global empty-interface case (nHat=0).
% % ADDED: signature checks + missing-field error checks.
% % ADDED: strict invalid-input rejection tests (char/logical/NaN/Inf/non-integer/duplicates/out-of-range).
%   If any strict rejection test fails with "DIAGNOSTIC: bad-case ... did not throw",
%   treat as Rule 7.3 and harden build_product_interface (do not relax the test).

  ensure_project_paths_();  % % CHANGED: now aligned with test_setup_local_schur.m style

  % ----------------------------
  % Case A: multi-subdomain array; one with empty glob_G while nHat>0
  % ----------------------------
  % Global DOFs: 1..6
  % Interface DOFs: 2,3,4,5  (hat indices 1..4)

  iface = struct();
  iface.nHat = 4;
  iface.glob2hat = zeros(6,1);
  iface.glob2hat(2) = 1;
  iface.glob2hat(3) = 2;
  iface.glob2hat(4) = 3;
  iface.glob2hat(5) = 4;

  sub = repmat(struct('glob_G', []), 5, 1);  % % ADDED: 5th has empty glob_G
  sub(1).glob_G = [2, 3];
  sub(2).glob_G = [3, 4];
  sub(3).glob_G = [2, 5];
  sub(4).glob_G = [4, 5];
  sub(5).glob_G = [];

  [sub_out, prod] = build_product_interface(sub, iface);

  % --- basic counts
  expected_nProd = 0;
  for i = 1:numel(sub)
    expected_nProd = expected_nProd + numel(sub(i).glob_G);
  end

  assert(prod.nHat  == iface.nHat,    sprintf('nHat mismatch: got %d expected %d', prod.nHat, iface.nHat));
  assert(prod.nProd == expected_nProd, sprintf('nProd mismatch: got %d expected %d', prod.nProd, expected_nProd));

  % sum_i |gamma_glob| == nProd
  sum_local = 0;
  for i = 1:numel(sub_out)
    sum_local = sum_local + numel(sub_out(i).gamma_glob);
  end
  assert(sum_local == prod.nProd, sprintf('Sum of local gamma sizes != nProd: %d vs %d', sum_local, prod.nProd));

  % --- product indices cover 1..nProd exactly once
  all_idx = [];
  for i = 1:numel(sub_out)
    all_idx = [all_idx; sub_out(i).prod_idx(:)]; %#ok<AGROW>
  end
  all_idx_sorted = sort(all_idx);
  assert(numel(all_idx_sorted) == prod.nProd, 'Total number of prod_idx entries != nProd');
  assert(all(all_idx_sorted == (1:prod.nProd).'), 'prod_idx are not a permutation of 1..nProd');

  % --- prod2sub/prod2hat ranges
  assert(all(prod.prod2sub >= 1 & prod.prod2sub <= numel(sub_out)), 'prod2sub contains out-of-range entries');
  assert(all(prod.prod2hat >= 1 & prod.prod2hat <= prod.nHat),        'prod2hat contains out-of-range entries');

  % --- per-subdomain consistency
  for i = 1:numel(sub_out)
    gi = sub_out(i).gamma_glob;
    hi = sub_out(i).gamma_hat;
    pi = sub_out(i).prod_idx;

    assert(numel(gi) == numel(hi), sprintf('sub(%d): gamma_glob and gamma_hat size mismatch', i));
    assert(numel(gi) == numel(pi), sprintf('sub(%d): gamma_glob and prod_idx size mismatch', i));

    if ~isempty(gi)
      mapped = iface.glob2hat(gi(:));
      assert(all(mapped == hi(:)), sprintf('sub(%d): gamma_hat != glob2hat(gamma_glob)', i));
      assert(all(mapped >= 1 & mapped <= iface.nHat), sprintf('sub(%d): gamma_hat contains invalid hat indices', i));
    end

    if ~isempty(pi)
      assert(all(prod.prod2sub(pi) == i), sprintf('sub(%d): prod2sub mismatch on prod_idx', i));
      assert(all(prod.prod2hat(pi) == hi(:)), sprintf('sub(%d): prod2hat mismatch on prod_idx', i));
    end
  end

  % % ADDED: explicit empty-subdomain check (nHat>0)
  assert(isempty(sub_out(5).gamma_glob) && isempty(sub_out(5).gamma_hat) && isempty(sub_out(5).prod_idx), ...
         'Expected empty gamma/prod_idx for subdomain with empty glob_G (nHat>0).');

  % --- multiplicity consistency between prod2hat and hat2prod
  mult_count = zeros(prod.nHat, 1);
  for j = 1:prod.nProd
    mult_count(prod.prod2hat(j)) = mult_count(prod.prod2hat(j)) + 1;
  end

  for hk = 1:prod.nHat
    list = prod.hat2prod{hk};
    assert(numel(list) == mult_count(hk), sprintf('hat2prod multiplicity mismatch for hk=%d', hk));
    if ~isempty(list)
      assert(all(prod.prod2hat(list) == hk), sprintf('hat2prod contains wrong indices for hk=%d', hk));
    end
  end

  % In this fixture each interface hat appears twice
  assert(all(mult_count == 2), sprintf('Expected all multiplicities 2, got [%s]', num2str(mult_count.')));

  % --- deterministic gather check
  u_hat  = [11; 22; 33; 44];
  u_prod = u_hat(prod.prod2hat);
  for hk = 1:prod.nHat
    list = prod.hat2prod{hk};
    for r = 1:numel(list)
      j = list(r);
      assert(u_prod(j) == u_hat(hk), sprintf('Gather mismatch: u_prod(%d) != u_hat(%d)', j, hk));
    end
  end

  % ----------------------------
  % Case B: global empty interface (nHat=0, all glob_G empty)
  % ----------------------------
  iface2 = struct('nHat', 0, 'glob2hat', zeros(0,1));
  sub2 = repmat(struct('glob_G', []), 3, 1);

  [sub2o, prod2] = build_product_interface(sub2, iface2);

  assert(prod2.nHat == 0 && prod2.nProd == 0, 'Expected (nHat,nProd)=(0,0) for empty case.');
  assert(isempty(prod2.prod2sub) && isempty(prod2.prod2hat), 'Expected empty prod2sub/prod2hat for empty case.');
  assert(iscell(prod2.hat2prod) && numel(prod2.hat2prod) == 0, 'Expected 0x1 cell hat2prod for empty case.');
  for i = 1:numel(sub2o)
    assert(isempty(sub2o(i).gamma_glob) && isempty(sub2o(i).gamma_hat) && isempty(sub2o(i).prod_idx), ...
      'Expected empty gamma/prod_idx for empty subdomains.');
  end

  % ----------------------------
  % Signature / error behavior
  % ----------------------------
  assert_throws_(@() build_product_interface(), 'build_product_interface nargin=0');                 % % ADDED
  assert_throws_(@() build_product_interface(sub, iface, 7), 'build_product_interface nargin=3');    % % ADDED
  assert_throws_(@() call_three_outputs_(sub, iface), 'build_product_interface too many outputs');   % % ADDED

  % ----------------------------
  % Missing required fields must error
  % ----------------------------
  assert_throws_(@() build_product_interface(sub, struct()), 'iface missing required fields');       % % ADDED

  tmp = iface; tmp = rmfield(tmp, 'glob2hat');
  assert_throws_(@() build_product_interface(sub, tmp), 'missing iface.glob2hat');                  % % ADDED

  tmp = iface; tmp = rmfield(tmp, 'nHat');
  assert_throws_(@() build_product_interface(sub, tmp), 'missing iface.nHat');                      % % ADDED

  tmpSub = repmat(struct(), 1, 1);
  assert_throws_(@() build_product_interface(tmpSub, iface), 'missing sub.glob_G');                 % % ADDED

  % ----------------------------
  % Strict invalid-input rejection (Rule 7.3)
  % ----------------------------
  % If any of these do NOT error, harden build_product_interface input validation
  % (reject char/logical; require numeric real finite; enforce index validity).

  bad = {};

  % % ADDED: char indexing hazard (make it in-range so it would "work" if accepted)
  iface_char = struct();
  iface_char.glob2hat = zeros(100,1);
  iface_char.glob2hat(65) = 1;   % 'A' -> 65 maps to hat 1
  iface_char.nHat = 1;

  s = struct('glob_G', 'A');
  bad{end+1} = struct('name', 'sub.glob_G is char (ASCII indexing hazard) [in-range sentinel]', ...
                      'sub', s, 'iface', iface_char);

  % % ADDED: logical mask hazard (size-compatible; would select valid entries if accepted)
  s = struct('glob_G', logical([0; 1; 0; 1; 0; 0])); % selects globals 2 and 4
  bad{end+1} = struct('name', 'sub.glob_G is logical mask (mask indexing hazard) [size-compatible]', ...
                      'sub', s, 'iface', iface);

  % non-integer
  s = struct('glob_G', [2; 3.5]);
  bad{end+1} = struct('name', 'sub.glob_G non-integer', 'sub', s, 'iface', iface);

  % out-of-range
  s = struct('glob_G', 999);
  bad{end+1} = struct('name', 'sub.glob_G out of range for iface.glob2hat', 'sub', s, 'iface', iface);

  % zero index
  s = struct('glob_G', 0);
  bad{end+1} = struct('name', 'sub.glob_G contains 0', 'sub', s, 'iface', iface);

  % duplicates
  s = struct('glob_G', [2; 2]);
  bad{end+1} = struct('name', 'sub.glob_G contains duplicates', 'sub', s, 'iface', iface);

  % non-interface dof (glob2hat==0)
  s = struct('glob_G', 1);
  bad{end+1} = struct('name', 'sub.glob_G includes non-interface dof (glob2hat=0)', 'sub', s, 'iface', iface);

  % NaN / Inf
  s = struct('glob_G', [2; NaN]);
  bad{end+1} = struct('name', 'sub.glob_G contains NaN', 'sub', s, 'iface', iface);

  s = struct('glob_G', [2; Inf]);
  bad{end+1} = struct('name', 'sub.glob_G contains Inf', 'sub', s, 'iface', iface);

  % invalid iface.glob2hat type
  iface_bad = iface; iface_bad.glob2hat = 'not-a-map';
  bad{end+1} = struct('name', 'iface.glob2hat is char', 'sub', sub(1), 'iface', iface_bad);

  for k = 1:numel(bad)
    threw = false;
    try
      build_product_interface(bad{k}.sub, bad{k}.iface);
    catch
      threw = true;
    end
    if ~threw
      fprintf('DIAGNOSTIC: bad-case %d did not throw (Rule 7.3): %s\n', k, bad{k}.name);
    end
    assert(threw, sprintf('Expected an error, but none was thrown (bad-case %d: %s).', k, bad{k}.name));
  end

  fprintf('PASS: test_build_product_interface\n');
end

% =========================
% Helpers
% =========================

function ensure_project_paths_()
% % ADDED: Make test runnable from any working directory (once tests/ is on path).
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

function assert_throws_(fh, label)
% Similar style to test_setup_local_schur.m: includes a diagnostic print.
  threw = false;
  try
    fh();
  catch
    threw = true;
  end

  if ~threw
    fprintf('DIAG: Expected an error but none was thrown for case: %s\n', label);
  end
  assert(threw, sprintf('Expected an error, but none was thrown (%s).', label));
end

function call_three_outputs_(sub, iface)
  [a, b, c] = build_product_interface(sub, iface); %#ok<NASGU,ASGLU>
end