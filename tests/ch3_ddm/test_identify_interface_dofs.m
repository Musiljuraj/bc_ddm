function test_identify_interface_dofs()
%TEST_IDENTIFY_INTERFACE_DOFS  Strong invariants for interior/interface split (Chapter 3.1.3).
%
% % CHANGED: the test now self-manages paths (can run from any folder),
% %          provided Octave can *find this test function* (see run instructions below).

  ensure_project_paths_();  % ADDED

  % ============================================================
  % Signature / contract checks
  % ============================================================

  % ADDED: wrong-arg-count must throw (do NOT match message substrings in Octave).
  assert_throws__(@() identify_interface_dofs(), ...
    'Expected an error for missing inputs.'); % ADDED
  assert_throws__(@() identify_interface_dofs(struct(), struct(), 1), ...
    'Expected an error for too many inputs.'); % ADDED

  % ADDED: too many outputs should throw (Octave-level error is OK).
  assert_throws__(@() eval('[a,b,c] = identify_interface_dofs(struct(), struct(''nFree'',1));'), ...
    'Expected an error for too many outputs.'); % ADDED

  % ============================================================
  % Minimal hand-constructed fixtures (cheap, deterministic)
  % ============================================================

  % ----------------------------
  % ADDED: One shared DOF across 3 subdomains
  % ----------------------------
  sub = struct('loc2glob', { [1;2;3], [3;4], [3;5] }); % ADDED
  ddm = struct('nFree', 5);                            % ADDED

  [sub2, iface] = identify_interface_dofs(sub, ddm);   % ADDED

  % ADDED: exact multiplicities
  assert(isequal(iface.counts(:), [1;1;3;1;1]), 'Bad counts in minimal shared-DOF fixture.'); % ADDED
  assert(isequal(iface.glob(:), 3), 'Bad iface.glob in minimal shared-DOF fixture.');       % ADDED
  assert(iface.nHat == 1, 'Bad iface.nHat in minimal shared-DOF fixture.');                 % ADDED

  % ADDED: mapping invariants
  assert(isequal(iface.hat2glob(:), iface.glob(:)), 'hat2glob must equal iface.glob.');     % ADDED
  assert(iface.glob2hat(3) == 1, 'glob2hat must map iface.glob -> 1:nHat.');                % ADDED
  assert(all(iface.glob2hat([1;2;4;5]) == 0), 'glob2hat must be zero off-interface.');      % ADDED
  assert(iface.glob2hat(iface.hat2glob(1)) == 1, 'glob2hat(hat2glob) must be identity.');  % ADDED

  % ADDED: dof2sub correctness
  assert(numel(iface.dof2sub) == iface.nHat, 'dof2sub size mismatch.');                     % ADDED
  assert(isequal(iface.dof2sub{1}(:), [1;2;3]), 'dof2sub wrong for minimal shared-DOF.');   % ADDED
  assert(numel(unique(iface.dof2sub{1})) == numel(iface.dof2sub{1}), 'dof2sub must not contain duplicates.'); % ADDED
  assert(numel(iface.dof2sub{1}) == iface.counts(3), 'dof2sub length must match multiplicity.'); % ADDED

  % ADDED: per-subdomain partition invariants
  for i = 1:numel(sub2)
    nLoc = numel(sub2(i).loc2glob(:));
    all_ids = sort([sub2(i).dofs_I(:); sub2(i).dofs_G(:)]);
    assert(isequal(all_ids, (1:nLoc).'), sprintf('Local dof partition not a full cover in sub %d.', i)); % ADDED
    assert(isempty(intersect(sub2(i).dofs_I(:), sub2(i).dofs_G(:))), ...
      sprintf('Local dof partition not disjoint in sub %d.', i)); % ADDED
    assert(sub2(i).nI == numel(sub2(i).dofs_I), sprintf('nI mismatch in sub %d.', i)); % ADDED
    assert(sub2(i).nG == numel(sub2(i).dofs_G), sprintf('nG mismatch in sub %d.', i)); % ADDED
  end

  % ----------------------------
  % ADDED: No shared DOFs edge case
  % ----------------------------
  sub = struct('loc2glob', { [1;2], [3;4] }); % ADDED
  ddm = struct('nFree', 4);                  % ADDED
  [sub2, iface] = identify_interface_dofs(sub, ddm); % ADDED

  assert(isequal(iface.counts(:), [1;1;1;1]), 'Counts must all be 1 when nothing is shared.'); % ADDED
  assert(isempty(iface.glob), 'iface.glob must be empty when nothing is shared.');             % ADDED
  assert(iface.nHat == 0, 'iface.nHat must be 0 when nothing is shared.');                     % ADDED
  assert(all(iface.glob2hat(:) == 0), 'glob2hat must be all zeros when no interface.');        % ADDED
  assert(isempty(iface.hat2glob), 'hat2glob must be empty when no interface.');                % ADDED
  assert(iscell(iface.dof2sub) && numel(iface.dof2sub) == 0, 'dof2sub must be empty cell when no interface.'); % ADDED

  for i = 1:numel(sub2)
    assert(isempty(sub2(i).dofs_G), sprintf('dofs_G must be empty in sub %d when no interface.', i)); % ADDED
    assert(isempty(sub2(i).glob_G), sprintf('glob_G must be empty in sub %d when no interface.', i)); % ADDED
    assert(sub2(i).nG == 0, sprintf('nG must be 0 in sub %d when no interface.', i)); % ADDED
  end

  % ============================================================
  % Small integration check (mesh + partition) but kept tiny
  % ============================================================

  n = 2;      % CHANGED
  nSubX = 2;  % CHANGED
  nSubY = 2;  % CHANGED

  [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);

  [sub, iface] = identify_interface_dofs(sub, ddm);

  % FIXED: verify multiplicity sum equals sum of local sizes (exact integer check).
  sum_local = 0;
  for i = 1:numel(sub)
    sum_local = sum_local + numel(sub(i).loc2glob);
  end
  assert(sum(iface.counts) == sum_local, 'sum(counts) must equal sum of local dof counts.'); % FIXED

  assert(all(iface.counts(iface.glob) >= 2), 'Interface DOFs must have multiplicity >= 2.');

  for i = 1:numel(sub)
    g = sub(i).loc2glob(:);
    maskG = iface.counts(g) > 1;
    assert(isequal(sub(i).glob_G(:), g(maskG)), sprintf('glob_G mismatch in sub %d.', i)); % CHANGED
    assert(isequal(sub(i).glob_I(:), g(~maskG)), sprintf('glob_I mismatch in sub %d.', i)); % CHANGED
  end

  assert(isequal(iface.hat2glob(:), iface.glob(:)), 'Integration: hat2glob must equal iface.glob.'); % ADDED
  if iface.nHat > 0
    assert(isequal(iface.glob2hat(iface.hat2glob(:)), (1:iface.nHat).'), ...
      'Integration: glob2hat(hat2glob) must equal 1:nHat.'); % ADDED
  end

  for k = 1:iface.nHat
    g = iface.hat2glob(k);
    expected = [];
    for i = 1:numel(sub)
      if any(sub(i).loc2glob(:) == g)
        expected(end+1,1) = i; %#ok<AGROW>
      end
    end
    got = iface.dof2sub{k}(:);
    assert(numel(unique(got)) == numel(got), sprintf('Integration: dof2sub{%d} has duplicates.', k)); % ADDED
    assert(isequal(sort(got), expected(:)), sprintf('Integration: dof2sub{%d} wrong membership.', k)); % ADDED
    assert(numel(got) == iface.counts(g), sprintf('Integration: dof2sub{%d} length must match multiplicity.', k)); % ADDED
  end

  % ============================================================
  % Invalid-input rejection (bug-finding). If these fail, harden the function.
  % ============================================================

  bad = {}; % ADDED
  bad{end+1} = struct('name','ddm missing nFree', ...
                      'sub', struct('loc2glob',{[1]}), ...
                      'ddm', struct()); % ADDED
  bad{end+1} = struct('name','loc2glob out of range (>nFree)', ...
                      'sub', struct('loc2glob',{[1;5]}), ...
                      'ddm', struct('nFree',4)); % ADDED
  bad{end+1} = struct('name','loc2glob is char (ASCII indices) should be rejected', ...
                      'sub', struct('loc2glob',{'a'}), ...
                      'ddm', struct('nFree',200)); % ADDED
  bad{end+1} = struct('name','loc2glob is logical should be rejected', ...
                      'sub', struct('loc2glob',{true}), ...
                      'ddm', struct('nFree',5)); % ADDED
  bad{end+1} = struct('name','ddm.nFree non-scalar should be rejected', ...
                      'sub', struct('loc2glob',{[1;2]}), ...
                      'ddm', struct('nFree',[4 1])); % ADDED

  for k = 1:numel(bad)
    did_throw = false;
    try
      identify_interface_dofs(bad{k}.sub, bad{k}.ddm);
    catch
      did_throw = true;
    end
    if ~did_throw
      error(sprintf('Expected an error but none was thrown for bad-case %d: %s', k, bad{k}.name)); % ADDED
    end
  end

  fprintf('PASS: test_identify_interface_dofs (tiny integration n=%d, nSubX=%d, nSubY=%d)\n', n, nSubX, nSubY);
end

% ============================================================
% Helpers
% ============================================================

function ensure_project_paths_()
% ADDED: locate project root from this test file and add main/tests paths

  % If setup_paths is already on path, prefer it (and still add tests tree).
  if exist('setup_paths','file') == 2
    sp = which('setup_paths');
    maindir = fileparts(sp);
    rootdir = fileparts(maindir);
    addpath(maindir);
    setup_paths();
    addpath(genpath(fullfile(rootdir, 'tests')));
    return;
  end

  % Otherwise: locate root relative to this file:
  %   .../tests/<chapter>/test_identify_interface_dofs.m
  thisdir  = fileparts(mfilename('fullpath'));
  testsdir = fileparts(thisdir);
  rootdir  = fileparts(testsdir);
  maindir  = fullfile(rootdir, 'main');

  if exist(fullfile(maindir, 'setup_paths.m'), 'file') ~= 2
    error('Could not locate main/setup_paths.m relative to test file location.');
  end

  addpath(maindir);
  setup_paths();

  % Add tests tree so helpers are found even when not run from root
  addpath(genpath(fullfile(rootdir, 'tests')));
end

function assert_throws__(fh, msg)
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  assert(threw, msg);
end