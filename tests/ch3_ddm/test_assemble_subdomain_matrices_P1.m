function test_assemble_subdomain_matrices_P1()
%TEST_ASSEMBLE_SUBDOMAIN_MATRICES_P1  Strengthened tests for local subdomain assembly (P1).
%
% Verifies that summing lifted local subdomain contributions reproduces the
% globally assembled Dirichlet-reduced system, without assuming identical DOF
% ordering between ddm and apply_dirichlet_elimination.
%
% Also checks basic invariants (sizes, sparsity, symmetry) and input contracts.

  local_setup_paths(); % ADDED (path-robust test runner)

  % Deterministic load (affine), should be handled robustly.
  f = @(x,y) 2.0 + x - 3*y;

  % Multiple deterministic cases (cheap meshes, divisibility respected).
  cases = { ...
    struct('n', 2, 'nSubX', 1, 'nSubY', 1), ...  % ADDED
    struct('n', 4, 'nSubX', 2, 'nSubY', 2), ...  % ADDED
    struct('n', 6, 'nSubX', 3, 'nSubY', 2)  ...  % CHANGED (keeps original coverage)
  };

  for ci = 1:numel(cases)
    c = cases{ci};
    run_consistency_case(c.n, c.nSubX, c.nSubY, f); % ADDED
  end

  fprintf('PASS: test_assemble_subdomain_matrices_P1 (%d cases)\n', numel(cases));
end


% -------------------------- helpers (local) -------------------------------

function run_consistency_case(n, nSubX, nSubY, f)
  % Build mesh and boundary info.
  [p, t, bnd] = mesh_unit_square_P1(n);

  % Global reference.
  K  = assemble_stiffness_P1(p, t);
  F  = assemble_load_P1(p, t, f);
  [Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, full(F), bnd.dirichlet_nodes);

  % Subdomain objects.
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
  [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<ASGLU>

  % --- Contract tests that need real objects available ---
  assert_throws(@() assemble_subdomain_matrices_P1(), 'nargin0'); % ADDED
  assert_throws(@() assemble_subdomain_matrices_P1(p), 'nargin1'); % ADDED
  assert_throws(@() assemble_subdomain_matrices_P1(p,t,sub), 'nargin3'); % ADDED
  assert_throws(@() assemble_subdomain_matrices_P1(p,t,sub,ddm,f,123), 'nargin6'); % ADDED

  assert_throws(@() assemble_subdomain_matrices_P1(p,t,sub,ddm,123), 'f_handle_not_function'); % ADDED

  ddm_bad = ddm;
  if isfield(ddm_bad, 'node2dof')
    ddm_bad = rmfield(ddm_bad, 'node2dof');
  end
  assert_throws(@() assemble_subdomain_matrices_P1(p,t,sub,ddm_bad,f), 'ddm_missing_node2dof'); % ADDED

  % Default f_handle should produce zero RHS.
  sub0 = assemble_subdomain_matrices_P1(p, t, sub, ddm); % ADDED
  for i = 1:numel(sub0)
    assert(issparse(sub0(i).f), 'Expected sub(i).f to be sparse.'); % ADDED
    assert(norm(full(sub0(i).f), 2) == 0, 'Default f_handle should yield zero RHS.'); % ADDED
  end

  % Assemble with nonzero f_handle.
  subF = assemble_subdomain_matrices_P1(p, t, sub, ddm, f);

  % Local invariants per subdomain.
  for i = 1:numel(subF)
    nloc = subF(i).nloc;

    assert(issparse(subF(i).K), 'Expected sub(i).K to be sparse.'); % ADDED
    assert(issparse(subF(i).f), 'Expected sub(i).f to be sparse.'); % ADDED

    assert(isequal(size(subF(i).K), [nloc, nloc]), ...
           sprintf('Size mismatch: sub(%d).K not nloc-by-nloc.', i)); % ADDED
    assert(isequal(size(subF(i).f), [nloc, 1]), ...
           sprintf('Size mismatch: sub(%d).f not nloc-by-1.', i)); % ADDED

    if nloc > 0
      % loc2glob sanity (should be a list of free DOFs).
      g = subF(i).loc2glob(:);
      assert(numel(g) == nloc, sprintf('sub(%d).loc2glob length != nloc.', i)); % ADDED
      assert(all(g > 0), sprintf('sub(%d).loc2glob contains non-positive entries.', i)); % ADDED
      assert(all(g == fix(g)), sprintf('sub(%d).loc2glob contains non-integers.', i)); % ADDED
      assert(all(g <= ddm.nFree), sprintf('sub(%d).loc2glob contains out-of-range DOFs.', i)); % ADDED

      % Symmetry of stiffness (Laplacian P1 stiffness should be symmetric).
      scaleKloc = max(1, norm(subF(i).K, 'fro'));
      tolSym = 2000 * eps * scaleKloc; % ADDED (scale-aware)
      assert(norm(subF(i).K - subF(i).K', 'fro') <= tolSym, ...
             sprintf('sub(%d).K is not symmetric within tolerance.', i)); % ADDED

      % Diagonal should be nonnegative for standard Poisson stiffness.
      d = diag(subF(i).K);
      assert(all(d >= -tolSym), sprintf('sub(%d).K has negative diagonal beyond tol.', i)); % ADDED
    end
  end

  % Lift + sum into global free-DOF space (in ddm DOF ordering).
  nFree = ddm.nFree;
  Ksum = sparse(nFree, nFree);
  Fsum = zeros(nFree, 1);

  for i = 1:ddm.Nsub
    g = subF(i).loc2glob(:);
    Ksum(g,g) = Ksum(g,g) + subF(i).K;
    Fsum(g)   = Fsum(g)   + full(subF(i).f);
  end

  % CHANGED: do NOT assume ddm.dof2node ordering equals free_nodes ordering.
  % Build a node->position map for the global reduced system, then permute Kff/Ff
  % into ddm's dof ordering.
  nNodes = size(p, 1);
  pos_in_free = zeros(nNodes, 1);                    % ADDED
  pos_in_free(free_nodes(:)) = (1:numel(free_nodes))'; % ADDED

  idx = pos_in_free(ddm.dof2node(:)); % ADDED
  assert(all(idx > 0), 'ddm.dof2node contains nodes not in free_nodes.'); % ADDED

  Kff_ddm = Kff(idx, idx); % ADDED
  Ff_ddm  = Ff(idx);       % ADDED

  % Scale-aware tolerances for global comparisons.
  scaleK = max(1, norm(Kff_ddm, 'fro')); % CHANGED
  scaleF = max(1, norm(Ff_ddm, 2));      % CHANGED
  tolK = 2000 * eps * scaleK;            % CHANGED
  tolF = 2000 * eps * scaleF;            % CHANGED

  assert(norm(Ksum - Kff_ddm, 'fro') <= tolK, ...
         sprintf('K mismatch (n=%d, %dx%d).', n, nSubX, nSubY)); % CHANGED
  assert(norm(Fsum - Ff_ddm, 2) <= tolF, ...
         sprintf('F mismatch (n=%d, %dx%d).', n, nSubX, nSubY)); % CHANGED

  % Order-independence sanity: reversing element order in one subdomain should not change
  % results beyond roundoff.
  subP = sub; % ADDED
  pick = 0;
  for i = 1:numel(subP)
    if subP(i).nloc > 0 && numel(subP(i).elems) > 1
      pick = i;
      break;
    end
  end
  if pick > 0
    subP(pick).elems = flipud(subP(pick).elems(:)); % ADDED
    subF2 = assemble_subdomain_matrices_P1(p, t, subP, ddm, f); % ADDED

    scaleKp = max(1, norm(subF(pick).K, 'fro'));
    scaleFp = max(1, norm(full(subF(pick).f), 2));
    tolKp = 5000 * eps * scaleKp; % ADDED
    tolFp = 5000 * eps * scaleFp; % ADDED

    assert(norm(subF2(pick).K - subF(pick).K, 'fro') <= tolKp, ...
           'Element-order change affected K beyond roundoff.'); % ADDED
    assert(norm(full(subF2(pick).f) - full(subF(pick).f), 2) <= tolFp, ...
           'Element-order change affected f beyond roundoff.'); % ADDED
  end

  % Contract enforcement: if a free DOF appears in elements but is removed from glob2loc, must error.
  subBad = sub; % ADDED
  bad_i = 0; bad_dof = 0; % ADDED
  for i = 1:numel(subBad)
    if subBad(i).nloc > 0 && numel(subBad(i).elems) > 0 && isnumeric(subBad(i).glob2loc)
      e = subBad(i).elems(1);
      nodes = t(e, :);
      dofs = ddm.node2dof(nodes);
      k = find(dofs > 0, 1);
      if ~isempty(k)
        bad_i = i;
        bad_dof = dofs(k);
        break;
      end
    end
  end
  if bad_i > 0 && bad_dof > 0
    if bad_dof <= numel(subBad(bad_i).glob2loc)
      subBad(bad_i).glob2loc(bad_dof) = 0; % ADDED
      assert_throws(@() assemble_subdomain_matrices_P1(p, t, subBad, ddm, f), ...
                    sprintf('glob2loc_missing_dof_i=%d', bad_i)); % ADDED
    end
  end

  % nloc==0 behavior: should produce empty sparse K and f, and not crash.
  extra = sub(1);      % ADDED
  extra.nloc = 0;      % ADDED
  extra.loc2glob = []; % ADDED
  extra.elems = [];    % ADDED

  subZ = sub(:);                 % FIXED (force column struct array)
  subZ(end+1,1) = extra;         % FIXED (append safely)
  subZ = assemble_subdomain_matrices_P1(p, t, subZ, ddm, f); % ADDED

  assert(isequal(size(subZ(end).K), [0 0]), 'nloc==0 should give 0x0 K.'); % ADDED
  assert(isequal(size(subZ(end).f), [0 1]), 'nloc==0 should give 0x1 f.'); % ADDED

  fprintf('  PASS case: n=%d, %dx%d subdomains\n', n, nSubX, nSubY); % CHANGED
end


function assert_throws(fcn, label)
  threw = false;
  try
    fcn();
  catch
    threw = true;
  end
  assert(threw, sprintf('Expected an error but none was thrown (%s).', label)); % ADDED
end


function local_setup_paths()
  thisfile = mfilename('fullpath');
  thisdir  = fileparts(thisfile);
  root = find_project_root(thisdir);

  addpath(fullfile(root, 'main'));                 % ADDED
  setup_paths();
  addpath(genpath(fullfile(root, 'tests')));       % ADDED
end


function root = find_project_root(startdir)
  root = startdir;
  for k = 1:10
    if exist(fullfile(root, 'main', 'setup_paths.m'), 'file')
      return;
    end
    parent = fileparts(root);
    if strcmp(parent, root)
      break;
    end
    root = parent;
  end
  error('Could not locate project root containing main/setup_paths.m'); % ADDED
end