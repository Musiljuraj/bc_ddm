function test_build_subdomains_structured()
%TEST_BUILD_SUBDOMAINS_STRUCTURED  Stronger checks for structured subdomain partitioning.
%
% % CHANGED: the test now self-manages paths (can run from any folder).

  % % FIXED: make setup_paths available even when running from tests/ch3_ddm
  ensure_project_paths_();  % ADDED

  % -----------------------------
  % 1) Positive cases (small meshes)
  % -----------------------------
  cases = {
    struct('n', 2, 'nSubX', 1, 'nSubY', 1)
    struct('n', 4, 'nSubX', 2, 'nSubY', 2)
    struct('n', 6, 'nSubX', 2, 'nSubY', 3)
  };

  for ci = 1:numel(cases)
    n     = cases{ci}.n;
    nSubX = cases{ci}.nSubX;
    nSubY = cases{ci}.nSubY;

    [p, t, bnd] = mesh_unit_square_P1(n);

    sub_only = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
    assert(isstruct(sub_only), sprintf('Expected struct output for sub (n=%d).', n));

    [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);

    assert(isstruct(ddm) && isscalar(ddm), sprintf('ddm must be a scalar struct (n=%d).', n));
    assert(ddm.nSubX == nSubX && ddm.nSubY == nSubY, sprintf('ddm nSub mismatch (n=%d).', n));
    assert(ddm.Nsub == nSubX*nSubY, sprintf('ddm.Nsub mismatch (n=%d).', n));

    nNodes = size(p,1);
    nElems = size(t,1);
    assert(ddm.nNodes == nNodes, sprintf('ddm.nNodes mismatch (n=%d).', n));
    assert(ddm.nElems == nElems, sprintf('ddm.nElems mismatch (n=%d).', n));

    dir_nodes = unique(bnd.dirichlet_nodes(:));
    all_nodes = (1:nNodes).';
    free_nodes = setdiff(all_nodes, dir_nodes);

    assert(isfield(ddm,'node2dof') && numel(ddm.node2dof) == nNodes, ...
           sprintf('ddm.node2dof size mismatch (n=%d).', n));
    assert(isfield(ddm,'dof2node') && numel(ddm.dof2node) == ddm.nFree, ...
           sprintf('ddm.dof2node size mismatch (n=%d).', n));
    assert(ddm.nFree == numel(free_nodes), sprintf('ddm.nFree mismatch (n=%d).', n));

    if ~isempty(dir_nodes)
      assert(all(ddm.node2dof(dir_nodes) == 0), sprintf('Dirichlet nodes must map to dof=0 (n=%d).', n));
    end
    if ~isempty(free_nodes)
      mapped = ddm.node2dof(free_nodes);
      assert(all(mapped >= 1 & mapped <= ddm.nFree), sprintf('Free nodes must map into 1..nFree (n=%d).', n));
      assert(numel(unique(mapped)) == numel(mapped), sprintf('node2dof must be injective on free nodes (n=%d).', n));
      assert(all(sort(mapped) == (1:ddm.nFree).'), sprintf('Free dofs must be exactly 1..nFree (n=%d).', n));
    end

    if ddm.nFree > 0
      invcheck = ddm.node2dof(ddm.dof2node(:));
      assert(all(invcheck == (1:ddm.nFree).'), sprintf('dof2node must invert node2dof (n=%d).', n));
    end

    assert(isstruct(sub) && numel(sub) == ddm.Nsub, sprintf('sub must have length Nsub (n=%d).', n));

    counts = zeros(nElems, 1);
    all_loc_dofs = [];

    nodes1 = t(:,1); nodes2 = t(:,2); nodes3 = t(:,3);
    cx = (p(nodes1,1) + p(nodes2,1) + p(nodes3,1)) / 3;
    cy = (p(nodes1,2) + p(nodes2,2) + p(nodes3,2)) / 3;
    tol = 50*eps;

    for s = 1:ddm.Nsub
      assert(sub(s).id == s, sprintf('sub(%d).id must equal %d (n=%d).', s, s, n));

      ix_s = mod(s-1, nSubX) + 1;
      iy_s = floor((s-1)/nSubX) + 1;
      assert(sub(s).ix == ix_s, sprintf('sub(%d).ix mismatch (n=%d).', s, n));
      assert(sub(s).iy == iy_s, sprintf('sub(%d).iy mismatch (n=%d).', s, n));

      x0 = (ix_s-1)/nSubX; x1 = ix_s/nSubX;
      y0 = (iy_s-1)/nSubY; y1 = iy_s/nSubY;

      assert(isfield(sub(s),'bbox') && numel(sub(s).bbox) == 4, sprintf('sub(%d).bbox malformed (n=%d).', s, n));
      bb = sub(s).bbox(:).';
      assert(all(bb == [x0 x1 y0 y1]), sprintf('sub(%d).bbox mismatch (n=%d).', s, n));

      elems = sub(s).elems(:);
      assert(~isempty(elems), sprintf('sub(%d) must have at least one element (n=%d).', s, n));
      assert(all(elems >= 1 & elems <= nElems), sprintf('sub(%d).elems out of range (n=%d).', s, n));
      counts(elems) = counts(elems) + 1;

      cxs = cx(elems);
      cys = cy(elems);
      assert(all(cxs >= x0 - tol & cxs <= x1 + tol), sprintf('sub(%d): centroid x outside bbox (n=%d).', s, n));
      assert(all(cys >= y0 - tol & cys <= y1 + tol), sprintf('sub(%d): centroid y outside bbox (n=%d).', s, n));

      nodes_expected = unique(t(elems,:));
      nodes_expected = nodes_expected(:);

      assert(isfield(sub(s),'nodes'), sprintf('sub(%d) missing nodes field (n=%d).', s, n));
      nodes_got = sub(s).nodes(:);
      assert(isequal(nodes_got, nodes_expected), sprintf('sub(%d).nodes inconsistent with its elements (n=%d).', s, n));

      dofs_expected = ddm.node2dof(nodes_expected);
      dofs_expected = dofs_expected(dofs_expected > 0);
      dofs_expected = unique(dofs_expected);
      dofs_expected = sort(dofs_expected);

      assert(isfield(sub(s),'loc2glob'), sprintf('sub(%d) missing loc2glob field (n=%d).', s, n));
      dofs_got = sub(s).loc2glob(:);
      assert(isequal(dofs_got, dofs_expected), sprintf('sub(%d).loc2glob mismatch vs expected (n=%d).', s, n));

      if ~isempty(dofs_got)
        assert(issorted(dofs_got), sprintf('sub(%d).loc2glob must be sorted (n=%d).', s, n));
        assert(all(dofs_got >= 1 & dofs_got <= ddm.nFree), sprintf('sub(%d).loc2glob out of range (n=%d).', s, n));
      end

      assert(isfield(sub(s),'nloc'), sprintf('sub(%d) missing nloc field (n=%d).', s, n));
      assert(sub(s).nloc == numel(dofs_got), sprintf('sub(%d).nloc mismatch (n=%d).', s, n));

      assert(isfield(sub(s),'glob2loc'), sprintf('sub(%d) missing glob2loc field (n=%d).', s, n));
      g2l = sub(s).glob2loc(:);
      assert(numel(g2l) == ddm.nFree, sprintf('sub(%d).glob2loc must be length nFree (n=%d).', s, n));

      if ~isempty(dofs_got)
        assert(all(g2l(dofs_got) == (1:numel(dofs_got)).'), sprintf('sub(%d).glob2loc inconsistent on support (n=%d).', s, n));
        assert(isequal(sub(s).loc2glob(g2l(dofs_got)), dofs_got), sprintf('sub(%d) loc2glob/glob2loc not inverse (n=%d).', s, n));
      end

      all_loc_dofs = [all_loc_dofs; dofs_got]; %#ok<AGROW>
    end

    assert(all(counts == 1), sprintf('Each element must belong to exactly one subdomain (n=%d).', n));

    if ddm.nFree > 0
      u = unique(all_loc_dofs(:));
      assert(isequal(u, (1:ddm.nFree).'), sprintf('Union of subdomain dofs must cover 1..nFree (n=%d).', n));
    end

    fprintf('PASS: test_build_subdomains_structured (n=%d, nSubX=%d, nSubY=%d)\n', n, nSubX, nSubY);
  end

  % -----------------------------
  % 2) Negative / contract tests
  % -----------------------------
  [p, t, bnd] = mesh_unit_square_P1(2);
  assert_throws(@() build_subdomains_structured(p, t, bnd, 1), ...
                'nargin=4 should throw');
  assert_throws(@() build_subdomains_structured(p, t, bnd, 1, 1, 7), ...
                'nargin=6 should throw');

  assert_throws(@() build_subdomains_structured(p(:,1), t, bnd, 1, 1), ...
                'p not Nx2 should throw');
  assert_throws(@() build_subdomains_structured(p, t(:,1:2), bnd, 1, 1), ...
                't not Nx3 should throw');

  bnd_missing = struct('dirichlet_nodes', bnd.dirichlet_nodes);
  assert_throws(@() build_subdomains_structured(p, t, bnd_missing, 1, 1), ...
                'bnd missing n should throw');
  bnd_missing = struct('n', bnd.n);
  assert_throws(@() build_subdomains_structured(p, t, bnd_missing, 1, 1), ...
                'bnd missing dirichlet_nodes should throw');

  assert_throws(@() build_subdomains_structured(p, t, bnd, 0, 1), ...
                'nSubX=0 should throw');
  assert_throws(@() build_subdomains_structured(p, t, bnd, 1, 0), ...
                'nSubY=0 should throw');
  assert_throws(@() build_subdomains_structured(p, t, bnd, 1.5, 1), ...
                'nSubX non-integer should throw');
  assert_throws(@() build_subdomains_structured(p, t, bnd, 1, -2), ...
                'nSubY negative should throw');

  bnd_bad = bnd;
  bnd_bad.n = 3;
  assert_throws(@() build_subdomains_structured(p, t, bnd_bad, 2, 1), ...
                'bnd.n not divisible by nSubX should throw');

  p_char = repmat('a', size(p,1), 2);
  assert_throws(@() build_subdomains_structured(p_char, t, bnd, 1, 1), ...
                'p as char matrix should throw (hardening target)');

  t_char = repmat('a', size(t,1), 3);
  assert_throws(@() build_subdomains_structured(p, t_char, bnd, 1, 1), ...
                't as char matrix should throw (hardening target)');

  bnd_oob = bnd;
  bnd_oob.dirichlet_nodes = [bnd.dirichlet_nodes(:); size(p,1)+5];
  assert_throws(@() build_subdomains_structured(p, t, bnd_oob, 1, 1), ...
                'out-of-range dirichlet_nodes should throw (hardening target)');

end

function ensure_project_paths_()
% % ADDED: locate project root from this test file and add main/tests paths
  if exist('setup_paths','file') == 2
    setup_paths();
    return;
  end

  thisdir = fileparts(mfilename('fullpath'));      % .../tests/ch3_ddm
  testsdir = fileparts(thisdir);                   % .../tests
  rootdir  = fileparts(testsdir);                  % .../ (project root)
  maindir  = fullfile(rootdir, 'main');

  if exist(fullfile(maindir, 'setup_paths.m'), 'file') ~= 2
    error('Could not locate main/setup_paths.m relative to test file location.');
  end

  addpath(maindir);
  setup_paths();

  % Add tests tree so helpers are found even when not run from root
  addpath(genpath(fullfile(rootdir, 'tests')));
end

function assert_throws(fh, tag)
  did_err = false;
  try
    fh();
  catch
    did_err = true;
  end
  assert(did_err, sprintf('Expected an error, but none was thrown: %s', tag));
end