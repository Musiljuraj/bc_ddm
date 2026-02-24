function ok = test_solve_tildeS()
%TEST_SOLVE_TILDES  Unit tests for solve_tildeS.m with explicit PASS/FAIL output
%
% Suggested location:
%   tests/ch4_feti_dp/test_solve_tildeS.m
%
% Run from project root:
%   addpath('main'); setup_paths(); addpath(genpath('tests')); test_solve_tildeS();

  ensure_project_paths_();

  test_name = mfilename();
  target_fn = 'solve_tildeS';

  fprintf('RUN : %s (target: %s)\n', test_name, target_fn);

  try
    ok = true;

    % ---------------------------
    % Case 1: small multi-subdomain synthetic data
    % ---------------------------
    data = struct();

    data.primal = struct('nC', 3);
    data.nDeltaProd = 3;

    sub = repmat(struct('c_ids', []), 3, 1);
    sub(1).c_ids = [1; 3];
    sub(2).c_ids = [2; 3];  % overlap at coarse id 3 (across subdomains)
    sub(3).c_ids = [1];     % primal-only subdomain
    data.sub = sub;

    data.delta_range = {1:2, 3:3, []};

    Sdd1 = [4, 1;
            1, 3];
    Sdd2 = 5;
    R1 = chol(Sdd1);  % upper, so R' * R = Sdd1
    R2 = chol(Sdd2);  % scalar
    data.Sdd_R = {R1, R2, []};

    Scd1 = [ 1.0,  0.5;
            -0.2,  0.3];     % (2 x 2)
    Scd2 = [ 0.7;
            -0.1];           % (2 x 1)
    data.Scd = {Scd1, Scd2, []};
    data.Sdc = {Scd1.', Scd2.', []};  % consistent synthetic choice

    Kcc = [6, 1, 0;
           1, 7, 1;
           0, 1, 5];
    data.Kcc_R = chol(Kcc);  % upper

    rc = [ 0.2; -1.0; 0.5];
    rD = [ 2.0; -3.0; 1.5];

    out = solve_tildeS(rc, rD, data);

    % Reference computation
    yD_ref = zeros(3,1);
    yD1 = R1 \ (R1' \ rD(1:2));
    yD2 = R2 \ (R2' \ rD(3));
    yD_ref(1:2) = yD1;
    yD_ref(3)   = yD2;

    rhs_c_ref = rc;
    rhs_c_ref([1;3]) = rhs_c_ref([1;3]) - Scd1 * yD1;
    rhs_c_ref([2;3]) = rhs_c_ref([2;3]) - Scd2 * yD2;

    Rcc = data.Kcc_R;
    u_c_ref = Rcc \ (Rcc' \ rhs_c_ref);

    wc_local_ref = cell(3,1);
    wc_local_ref{1} = u_c_ref([1;3]);
    wc_local_ref{2} = u_c_ref([2;3]);
    wc_local_ref{3} = u_c_ref([1]);

    wD_ref = yD_ref;
    wD_ref(1:2) = wD_ref(1:2) - (R1 \ (R1' \ (data.Sdc{1} * wc_local_ref{1})));
    wD_ref(3)   = wD_ref(3)   - (R2 \ (R2' \ (data.Sdc{2} * wc_local_ref{2})));

    assert(isstruct(out), 'Expected output to be a struct.');
    assert(isfield(out,'u_c') && isfield(out,'wc_local') && isfield(out,'wD'), ...
           'Output struct missing expected fields.');

    assert_close_(out.u_c, u_c_ref, tol_for_(u_c_ref), 'u_c mismatch.');
    assert_close_(out.wD,  wD_ref,  tol_for_(wD_ref),  'wD mismatch.');

    assert(iscell(out.wc_local) && numel(out.wc_local) == 3, 'wc_local should be a 3-cell.');
    for i = 1:3
      assert_close_(out.wc_local{i}, wc_local_ref{i}, tol_for_(wc_local_ref{i}), ...
                    sprintf('wc_local{%d} mismatch.', i));
    end

    % Defaulting behavior: empty rc/rD -> zeros of correct size
    out0 = solve_tildeS([], [], data);
    assert_close_(out0.u_c, zeros(3,1), tol_for_(zeros(3,1)), 'Default u_c not zero.');
    assert_close_(out0.wD,  zeros(3,1), tol_for_(zeros(3,1)), 'Default wD not zero.');
    for i = 1:3
      assert_close_(out0.wc_local{i}, zeros(numel(data.sub(i).c_ids),1), ...
                    tol_for_(zeros(numel(data.sub(i).c_ids),1)), ...
                    sprintf('Default wc_local{%d} not zero.', i));
    end

    % Linearity check
    rc1 = [1;0;-1];   rD1 = [0.5; 1.0; -0.5];
    rc2 = [-2;3; 0];  rD2 = [1.5; -2.0; 2.5];

    o1 = solve_tildeS(rc1, rD1, data);
    o2 = solve_tildeS(rc2, rD2, data);
    os = solve_tildeS(rc1+rc2, rD1+rD2, data);

    assert_close_(os.u_c, o1.u_c + o2.u_c, tol_for_(os.u_c), 'Linearity failed for u_c.');
    assert_close_(os.wD,  o1.wD  + o2.wD,  tol_for_(os.wD),  'Linearity failed for wD.');
    for i = 1:3
      assert_close_(os.wc_local{i}, o1.wc_local{i} + o2.wc_local{i}, tol_for_(os.wc_local{i}), ...
                    sprintf('Linearity failed for wc_local{%d}.', i));
    end

    % ---------------------------
    % Case 2: nC == 0 (delta-only)
    % ---------------------------
    data2 = struct();
    data2.primal = struct('nC', 0);
    data2.nDeltaProd = 1;
    data2.sub = struct('c_ids', {[]});
    data2.delta_range = {1:1};

    Sdd = 2;
    data2.Sdd_R = {chol(Sdd)};
    data2.Scd = {zeros(0,1)};
    data2.Sdc = {zeros(1,0)};
    data2.Kcc_R = [];

    out2 = solve_tildeS([], 1, data2);
    wD2_ref = (data2.Sdd_R{1} \ (data2.Sdd_R{1}' \ 1));
    assert(isempty(out2.u_c), 'Expected empty u_c when nC==0.');
    assert_close_(out2.wD, wD2_ref, tol_for_(wD2_ref), 'nC==0 wD mismatch.');

    % ---------------------------
    % Negative tests: wrong arg count + inconsistent indexing should throw
    % ---------------------------
    assert_throws_(@() solve_tildeS(), 'Expected error for 0 args.');
    assert_throws_(@() solve_tildeS(rc), 'Expected error for 1 arg.');
    assert_throws_(@() solve_tildeS(rc, rD), 'Expected error for 2 args.');

    data_bad = data;
    data_bad.delta_range{1} = 1:5; % out of bounds for rD length 3
    assert_throws_(@() solve_tildeS(rc, rD, data_bad), 'Expected error for out-of-range delta_range.');

    data_bad2 = data;
    data_bad2.sub(1).c_ids = [1; 4]; % out of bounds for nC==3
    assert_throws_(@() solve_tildeS(rc, rD, data_bad2), 'Expected error for out-of-range c_ids.');

    fprintf('PASS: %s (tested: %s)\n', test_name, target_fn);

  catch err
    ok = false;
    fprintf(2, 'FAIL: %s (tested: %s)\n', test_name, target_fn);
    fprintf(2, '  %s\n', err.message);
    rethrow(err);
  end
end


% ===== helper: tolerance =====
function t = tol_for_(x)
  nx = norm(x(:), 2);
  t = 200 * eps(max(1, nx));
end

% ===== helper: numeric closeness =====
function assert_close_(a, b, tol, msg)
  a = a(:); b = b(:);
  if numel(a) ~= numel(b)
    error('assert_close_: size mismatch. %s', msg);
  end
  d = norm(a - b, 2);
  if d > tol
    error('assert_close_: %s (||a-b||=%g > tol=%g)', msg, d, tol);
  end
end

% ===== helper: assert error thrown =====
function assert_throws_(f, msg)
  threw = false;
  try
    f();
  catch
    threw = true;
  end
  if ~threw
    error('assert_throws_: %s', msg);
  end
end

% ===== mandatory helper: ensure project paths =====
function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did
    return;
  end

  % Strategy 1: setup_paths already resolvable
  if exist('setup_paths','file') == 2
    setup_paths();
    sp = which('setup_paths');
    main_dir = fileparts(sp);
    root_dir = fileparts(main_dir);
    addpath(main_dir);
    addpath(genpath(fullfile(root_dir,'tests')));
    did = true;
    return;
  end

  % Strategy 2: derive root from this test file location, walk upwards
  thisdir = fileparts(mfilename('fullpath'));
  root_dir = thisdir;

  for k = 1:10
    cand = fullfile(root_dir, 'main', 'setup_paths.m');
    if exist(cand, 'file') == 2
      addpath(fullfile(root_dir, 'main'));
      setup_paths();
      addpath(genpath(fullfile(root_dir,'tests')));
      did = true;
      return;
    end
    parent = fileparts(root_dir);
    if strcmp(parent, root_dir)
      break;
    end
    root_dir = parent;
  end

  error('ensure_project_paths_: could not locate project root containing main/setup_paths.m');
end