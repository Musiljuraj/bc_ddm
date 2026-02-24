function test_select_primal_dofs()
%TEST_SELECT_PRIMAL_DOFS  Corner-based primal selection: invariants + input contract.
%
% % CHANGED: Path setup now uses ensure_project_paths_() (same pattern as test_setup_local_schur.m).
% % CHANGED: Expected primal set computed independently via enumerated corner coordinates.
% % ADDED  : Signature checks (wrong-arg-count, too-many-outputs).
% % ADDED  : Ordering invariance under permutation of iface.glob (with consistent glob2hat rebuild).
% % ADDED  : Input-rejection tests (likely to FAIL until select_primal_dofs is hardened).

  ensure_project_paths_();  % % CHANGED

  % Keep partition lines on mesh nodes: choose n divisible by nSubX and nSubY.
  cases = {                                                  % CHANGED
    struct('n',  4, 'nSubX', 2, 'nSubY', 2)
    struct('n',  6, 'nSubX', 3, 'nSubY', 2)
    struct('n', 12, 'nSubX', 4, 'nSubY', 3)
  };

  tol = 1e-12;  % should match current select_primal_dofs implementation

  % ---------- Run correctness + invariants across cases ----------
  for ci = 1:numel(cases)
    n     = cases{ci}.n;
    nSubX = cases{ci}.nSubX;
    nSubY = cases{ci}.nSubY;

    [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
    [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY); %#ok<ASGLU>
    [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>

    primal = select_primal_dofs(p, ddm, iface);

    % ---- Output shape/type invariants ----
    assert(isstruct(primal), 'primal must be a struct.');                         % ADDED
    local_assert_has_fields(primal, {'glob_c','glob_delta','hat_c','hat_delta','nC'}); % ADDED

    local_assert_int_vector(primal.glob_c,     'primal.glob_c');                  % ADDED
    local_assert_int_vector(primal.glob_delta, 'primal.glob_delta');              % ADDED
    local_assert_int_vector(primal.hat_c,      'primal.hat_c');                   % ADDED
    local_assert_int_vector(primal.hat_delta,  'primal.hat_delta');               % ADDED
    assert(isscalar(primal.nC) && primal.nC == numel(primal.glob_c), ...
           'primal.nC must equal numel(primal.glob_c).');                         % ADDED

    assert(isequal(primal.glob_c(:),     sort(primal.glob_c(:))),     'glob_c must be sorted.');     % ADDED
    assert(isequal(primal.glob_delta(:), sort(primal.glob_delta(:))), 'glob_delta must be sorted.'); % ADDED

    iface_glob = iface.glob(:);
    assert(numel(unique(iface_glob)) == numel(iface_glob), 'iface.glob must contain unique DOFs.');  % ADDED

    % ---- Partition invariants ----
    assert(isempty(intersect(primal.glob_c, primal.glob_delta)), 'glob_c and glob_delta must be disjoint.');
    union_sorted = sort([primal.glob_c(:); primal.glob_delta(:)]);
    assert(isequal(union_sorted, sort(iface_glob)), 'glob_c U glob_delta must equal iface.glob.');

    % ---- Hat mapping invariants ----
    nHat = numel(iface_glob);
    assert(all(primal.hat_c     >= 1 & primal.hat_c     <= nHat), 'hat_c out of range.');            % ADDED
    assert(all(primal.hat_delta >= 1 & primal.hat_delta <= nHat), 'hat_delta out of range.');        % ADDED
    assert(isempty(intersect(primal.hat_c, primal.hat_delta)), 'hat_c and hat_delta must be disjoint.'); % ADDED

    % Consistency: hat indices map back to the selected globals.
    assert(isequal(iface.glob(primal.hat_c(:)), primal.glob_c(:)), ...
           'iface.glob(hat_c) must equal glob_c (in the same order).');                               % ADDED
    assert(isequal(iface.glob(primal.hat_delta(:)), primal.glob_delta(:)), ...
           'iface.glob(hat_delta) must equal glob_delta (in the same order).');                       % ADDED

    % Also check direct glob2hat usage for the chosen sets (function contract).
    assert(isequal(primal.hat_c(:),     iface.glob2hat(primal.glob_c(:))),     'hat_c must equal glob2hat(glob_c).');
    assert(isequal(primal.hat_delta(:), iface.glob2hat(primal.glob_delta(:))), 'hat_delta must equal glob2hat(glob_delta).');

    % ---- Independent expected primal set via enumerated corner coordinates ----
    expected_glob_c = local_expected_primal_from_corners(p, ddm, iface, tol);     % CHANGED
    expected_glob_delta = sort(setdiff(iface_glob, expected_glob_c));             % ADDED

    assert(isequal(primal.glob_c(:), expected_glob_c(:)), 'glob_c mismatch vs enumerated corners.');      % CHANGED
    assert(isequal(primal.glob_delta(:), expected_glob_delta(:)), 'glob_delta mismatch vs complement.');  % ADDED

    % ---- Ordering invariance under permutation of iface.glob ----
    perm = (numel(iface_glob):-1:1)';                                             % ADDED (deterministic)
    iface2 = struct();
    iface2.glob = iface_glob(perm);                                               % ADDED

    % Rebuild a consistent glob2hat for the permuted ordering.
    iface2.glob2hat = zeros(size(iface.glob2hat));                                % ADDED
    for j = 1:numel(iface2.glob)
      g = iface2.glob(j);
      iface2.glob2hat(g) = j;
    end

    primal2 = select_primal_dofs(p, ddm, iface2);                                 % ADDED
    assert(isequal(primal2.glob_c(:), primal.glob_c(:)), 'glob_c must be invariant to iface ordering.');         % ADDED
    assert(isequal(primal2.glob_delta(:), primal.glob_delta(:)), 'glob_delta must be invariant to iface ordering.'); % ADDED
    assert(isequal(iface2.glob(primal2.hat_c(:)), primal2.glob_c(:)), 'Permuted hat_c must map back correctly.');     % ADDED
    assert(isequal(iface2.glob(primal2.hat_delta(:)), primal2.glob_delta(:)), 'Permuted hat_delta must map back correctly.'); % ADDED

    fprintf('PASS: test_select_primal_dofs (n=%d, %dx%d subdomains, nC=%d)\n', ...
            n, nSubX, nSubY, primal.nC);
  end

  % ---------- Signature checks (use a valid fixture) ----------
  n = 4; nSubX = 2; nSubY = 2;
  [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY); %#ok<ASGLU>
  [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>

  local_assert_throws_any(@() select_primal_dofs(),                    'nargin=0');  % ADDED
  local_assert_throws_any(@() select_primal_dofs(p),                   'nargin=1');  % ADDED
  local_assert_throws_any(@() select_primal_dofs(p, ddm),              'nargin=2');  % ADDED
  local_assert_throws_any(@() select_primal_dofs(p, ddm, iface, 123),  'nargin=4');  % ADDED

  % Too many outputs should throw.
  local_assert_throws_any(@() eval('[a,b] = select_primal_dofs(p, ddm, iface);'), 'nargout=2'); % ADDED

  % ---------- Input rejection checks (expected to FAIL until hardening) ----------
  primal_ok = select_primal_dofs(p, ddm, iface); %#ok<NASGU>

  bad_cases = {
    struct('label','p_char_matrix', ...
           'fh', @() select_primal_dofs(repmat('A', size(p)), ddm, iface))                     % ADDED
    struct('label','p_has_NaN', ...
           'fh', @() select_primal_dofs(local_set_one_entry(p, NaN), ddm, iface))             % ADDED
    struct('label','p_wrong_cols_3', ...
           'fh', @() select_primal_dofs([p, zeros(size(p,1),1)], ddm, iface))                 % ADDED
    struct('label','ddm_nSubX_zero', ...
           'fh', @() select_primal_dofs(p, local_setfield(ddm, 'nSubX', 0), iface))           % ADDED
    struct('label','ddm_nSubX_nonint', ...
           'fh', @() select_primal_dofs(p, local_setfield(ddm, 'nSubX', 2.5), iface))         % ADDED
    struct('label','iface_glob2hat_all_zero', ...
           'fh', @() select_primal_dofs(p, ddm, local_setfield(iface, 'glob2hat', zeros(size(iface.glob2hat))))) % ADDED
    struct('label','iface_glob_out_of_range', ...
           'fh', @() select_primal_dofs(p, ddm, local_setfield(iface, 'glob', [iface.glob(:); numel(ddm.dof2node)+1]))) % ADDED
    struct('label','ddm_missing_field', ...
           'fh', @() select_primal_dofs(p, rmfield(ddm, 'dof2node'), iface))                  % ADDED (should already throw)
    struct('label','iface_missing_field', ...
           'fh', @() select_primal_dofs(p, ddm, rmfield(iface, 'glob2hat')))                  % ADDED (should already throw)
  };

  for i = 1:numel(bad_cases)
    local_assert_throws_any(bad_cases{i}.fh, bad_cases{i}.label);                               % ADDED
  end

  fprintf('PASS: test_select_primal_dofs (signature + input-contract checks)\n');
end


% ======================= helpers =======================

function ensure_project_paths_()
% % CHANGED: Same path-robust helper pattern as test_setup_local_schur.m.
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

function local_assert_has_fields(s, fields)
% % ADDED
  for k = 1:numel(fields)
    f = fields{k};
    assert(isfield(s, f), sprintf('Missing required field: %s', f));
  end
end

function local_assert_int_vector(v, name)
% % ADDED: integer-valued, real, finite vector (empty allowed).
  assert(isvector(v), sprintf('%s must be a vector.', name));
  if isempty(v)
    return;
  end
  assert(isnumeric(v) && isreal(v), sprintf('%s must be real numeric.', name));
  assert(all(isfinite(v(:))), sprintf('%s must be finite.', name));
  assert(all(v(:) == round(v(:))), sprintf('%s must be integer-valued.', name));
end

function local_assert_throws_any(fh, label)
% % ADDED: Diagnostic for "expected error but none thrown".
  did_throw = false;
  try
    fh();
  catch
    did_throw = true;
  end
  if ~did_throw
    fprintf(2, 'NO ERROR THROWN for bad case: %s\n', label); % diagnostic
  end
  assert(did_throw, sprintf('Expected an error but none was thrown (%s).', label));
end

function expected_glob_c = local_expected_primal_from_corners(p, ddm, iface, tol)
% % CHANGED: Enumerate all decomposition corner coordinates and map back to free DOFs.
  nSubX = ddm.nSubX;
  nSubY = ddm.nSubY;

  xLines = (0:nSubX) / nSubX;
  yLines = (0:nSubY) / nSubY;

  iface_glob = iface.glob(:);
  dof2node = ddm.dof2node(:);

  g_list = [];

  for ix = 1:numel(xLines)
    for iy = 1:numel(yLines)
      x = xLines(ix);
      y = yLines(iy);

      node = find(abs(p(:,1) - x) < tol & abs(p(:,2) - y) < tol);
      assert(numel(node) == 1, sprintf('Expected exactly 1 node at (%.17g, %.17g).', x, y));

      dof = find(dof2node == node);
      assert(numel(dof) <= 1, 'ddm.dof2node must map each node to at most one free DOF.');

      if ~isempty(dof)
        g = dof(1);
        if ismember(g, iface_glob)
          g_list(end+1,1) = g; %#ok<AGROW>
        end
      end
    end
  end

  expected_glob_c = unique(sort(g_list));
end

function p2 = local_set_one_entry(p, val)
% % ADDED
  p2 = p;
  if ~isempty(p2)
    p2(1,1) = val;
  end
end

function s2 = local_setfield(s, name, value)
% % ADDED
  s2 = s;
  s2.(name) = value;
end