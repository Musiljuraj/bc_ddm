function test_build_problem_data()
%TEST_BUILD_PROBLEM_DATA  Contract/invariant tests for build_problem_data().
%
% Suggested location:
%   tests/ch4_common/test_build_problem_data.m

  ensure_project_paths_();

  try
    run_core_();
    fprintf('test_build_problem_data: PASS\n');
  catch ME
    if isfield(ME, 'identifier') && ~isempty(ME.identifier)
      fprintf('test_build_problem_data: FAIL (%s)\n', ME.identifier);
    else
      fprintf('test_build_problem_data: FAIL\n');
    end
    rethrow(ME);
  end
end

% ========================= core test =========================

function run_core_()
  % --- Smoke + structural invariants (tiny but VALID DDM case: >=2 subdomains) ---
  n = 2;
  nSubX = 2;
  nSubY = 1;

  data = build_problem_data(n, nSubX, nSubY, []);  % default f_handle path

  assert(isstruct(data));
  assert_struct_has_fields_(data, {'p','t','bnd','K','F','Kff','Ff','free','ddm','iface','sub','prod','primal'});

  % Mesh basics
  p = data.p; t = data.t;
  assert(isnumeric(p) && isreal(p) && size(p,2) == 2);
  assert(isnumeric(t) && isreal(t) && size(t,2) == 3);
  nNodes = size(p,1);
  assert(nNodes == (n+1)*(n+1));

  % Unit square coordinates
  assert(all(p(:) >= -1e-14) && all(p(:) <= 1+1e-14));

  % Connectivity indices
  assert(all(abs(t(:) - round(t(:))) < 1e-12));
  assert(all(t(:) >= 1) && all(t(:) <= nNodes));

  % Global K,F shapes
  K = data.K; F = data.F;
  assert(ismatrix(K) && size(K,1) == nNodes && size(K,2) == nNodes);
  assert(isvector(F) && numel(F) == nNodes);

  % Dirichlet elimination invariants
  assert(isfield(data.bnd, 'dirichlet_nodes'));
  dirn = data.bnd.dirichlet_nodes(:);
  assert(isnumeric(dirn) && isreal(dirn));
  assert(all(abs(dirn - round(dirn)) < 1e-12));
  assert(all(dirn >= 1) && all(dirn <= nNodes));
  dirn = unique(dirn);

  free = data.free(:);
  assert(isnumeric(free) && isreal(free));
  assert(all(abs(free - round(free)) < 1e-12));
  assert(all(free >= 1) && all(free <= nNodes));
  free = unique(free);

  comp = setdiff((1:nNodes).', dirn);
  assert(isequal(free, comp));

  Kff = data.Kff; Ff = data.Ff;
  assert(ismatrix(Kff) && size(Kff,1) == numel(free) && size(Kff,2) == numel(free));
  assert(isvector(Ff) && numel(Ff) == numel(free));

  % Recompute elimination from returned (K,F,bnd) and compare
  [Kff2, Ff2, free2] = apply_dirichlet_elimination(K, F, dirn);
  assert(isequal(free2(:), free));
  assert_close_(Kff, Kff2, 1e-12);
  assert_close_(Ff(:), Ff2(:), 1e-12);

  % Symmetry + SPD check on Kff
  sym_err = norm(Kff - Kff.', 'fro');
  sym_tol = 1e-12 * max(1, norm(Kff, 'fro'));
  assert(sym_err <= sym_tol);
  A = 0.5 * (Kff + Kff.');
  [~, pchol] = chol(A);
  assert(pchol == 0);

  % Interface should exist for >=2 subdomains
  assert(isfield(data.iface, 'glob'));
  assert(isvector(data.iface.glob) && ~isempty(data.iface.glob));

  % --- Default RHS equivalence: empty f_handle should match f(x,y)=1 ---
  data1 = build_problem_data(2, 2, 1, []);
  data2 = build_problem_data(2, 2, 1, @(x,y) 1.0);
  assert_close_(data1.F(:),   data2.F(:),   1e-13);
  assert_close_(data1.Ff(:),  data2.Ff(:),  1e-13);
  assert_close_(data1.K(:),   data2.K(:),   1e-13);
  assert_close_(data1.Kff(:), data2.Kff(:), 1e-13);

  % --- Partition count sanity (structured partition) ---
  dataP = build_problem_data(4, 2, 2, []);
  assert(isstruct(dataP.sub));
  assert(numel(dataP.sub) == 4);

  % --- Negative tests: wrong arg count ---
  assert_throws_(@() build_problem_data(), 'nargin0');
  assert_throws_(@() build_problem_data(2), 'nargin1');
  assert_throws_(@() build_problem_data(2,1), 'nargin2');
  assert_throws_(@() build_problem_data(2,1,1,[],5), 'nargin5');

  % --- Negative tests: 1x1 decomposition should throw (empty interface) ---
  assert_throws_(@() build_problem_data(2, 1, 1, []), 'singleSubdomain');

  % --- Negative tests: invalid types / values (should throw) ---
  bad = {
    @() build_problem_data('2', 1, 1, []),   'char_n';
    @() build_problem_data(NaN, 1, 1, []),   'nan_n';
    @() build_problem_data(2, 0, 1, []),     'zero_nSubX';
    @() build_problem_data(2, 1, -1, []),    'neg_nSubY';
    @() build_problem_data(3, 2, 1, []),     'nondividing_partition';
    @() build_problem_data(2, 1, 1, 7),      'nonhandle_rhs';
  };
  assert_throws_cases_(bad);
end

% ========================= helpers =========================

function assert_struct_has_fields_(s, fields)
  for k = 1:numel(fields)
    assert(isfield(s, fields{k}), ['missing field: ', fields{k}]);
  end
end

function assert_close_(A, B, reltol)
  if nargin < 3, reltol = 1e-12; end
  A = A(:); B = B(:);
  assert(numel(A) == numel(B));
  scale = max(1, max(norm(A), norm(B)));
  err = norm(A - B);
  assert(err <= reltol * scale);
end

function assert_throws_(fh, label)
  didThrow = false;
  try
    fh();
  catch
    didThrow = true;
  end
  if ~didThrow
    fprintf('Expected an error but none was thrown (%s)\n', label);
  end
  assert(didThrow);
end

function assert_throws_cases_(cases)
  for i = 1:size(cases,1)
    fh = cases{i,1};
    label = cases{i,2};
    assert_throws_(fh, label);
  end
end

function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did, return; end

  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    main_dir = fileparts(sp);
    root_dir = fileparts(main_dir);
    addpath(genpath(fullfile(root_dir, 'tests')));
    did = true;
    return;
  end

  oldpwd = pwd();
  c = onCleanup(@() cd(oldpwd)); %#ok<NASGU>

  thisdir = fileparts(mfilename('fullpath'));
  root_dir = thisdir;

  found = false;
  for k = 1:10
    cand = root_dir;
    if exist(fullfile(cand, 'main', 'setup_paths.m'), 'file') == 2
      root_dir = cand;
      found = true;
      break;
    end
    parent = fileparts(cand);
    if strcmp(parent, cand), break; end
    root_dir = parent;
  end

  if ~found
    error('ensure_project_paths_:notFound', ...
          'could not locate project root containing main/setup_paths.m');
  end

  addpath(fullfile(root_dir, 'main'));
  setup_paths();
  addpath(genpath(fullfile(root_dir, 'tests')));
  did = true;
end