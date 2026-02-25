function test_reconstruct_bddc_solution()
% Minimal-standard test for reconstruct_bddc_solution()

  fprintf('Running test_reconstruct_bddc_solution (function: reconstruct_bddc_solution)\n');

  % --- Noninteractive cleanup: restore path + clear cached functions ---
  orig_path = path();
  c_path = onCleanup(@() cleanup_paths_and_cache_(orig_path)); %#ok<NASGU>

  ensure_project_paths_();

  % --- wrong-arg-count must throw ---
  assert_throws_(@() reconstruct_bddc_solution(), 'no-args');
  assert_throws_(@() reconstruct_bddc_solution(1), 'one-arg');
  assert_throws_(@() reconstruct_bddc_solution(1, struct(), 3), 'three-args');

  % --- representative invalid input must throw (missing required field) ---
  data = build_tiny_data_();
  data2 = rmfield(data, 'bddc');
  assert_throws_(@() reconstruct_bddc_solution([1;2], data2), 'missing-field-bddc');

  % --- representative invalid input must throw (dimension mismatch in R*u_hat) ---
  assert_throws_(@() reconstruct_bddc_solution([1;2;3], data), 'u_hat-length-mismatch');

  % --- representative invalid input must throw (zero contributions to some free dof) ---
  data3 = data;
  data3.ddm.nFree = 3; % dof #3 never receives contributions => must throw
  assert_throws_(@() reconstruct_bddc_solution([2;3], data3), 'zero-contrib-free-dof');

  % --- smoke path + output contract checks ---
  u_hat = [2; 3];
  out = reconstruct_bddc_solution(u_hat, data);

  assert(isstruct(out));
  assert(isfield(out, 'w_prod'));
  assert(isfield(out, 'u_free'));
  assert(isfield(out, 'count'));

  % Shapes / sizes
  assert(isvector(out.w_prod));
  assert(numel(out.w_prod) == size(data.bddc.R, 1));

  assert(isvector(out.u_free));
  assert(numel(out.u_free) == data.ddm.nFree);

  assert(isvector(out.count));
  assert(numel(out.count) == data.ddm.nFree);

  % Real/finite numeric values where expected
  assert(isreal(out.w_prod) && all(isfinite(out.w_prod)));
  assert(isreal(out.u_free) && all(isfinite(out.u_free)));
  assert(isreal(out.count) && all(isfinite(out.count)));

  % Stable invariant: w_prod equals R*u_hat
  w_ref = data.bddc.R * u_hat(:);
  assert(norm(out.w_prod - w_ref, 2) <= 1e-12 * max(1, norm(w_ref, 2)));

  % Core sanity checks (tiny deterministic reference values)
  % - u_free(1) is average of the two product interface contributions (2 and 3) => 2.5
  % - u_free(2) is interior dof from sub(1): uI = (f_I - K_Ig*w_g) = 10 - 1*2 = 8
  assert(abs(out.u_free(1) - 2.5) <= 1e-12);
  assert(abs(out.u_free(2) - 8.0) <= 1e-12);

  % Multiplicity count check: dof #1 hit twice, dof #2 hit once
  assert(abs(out.count(1) - 2) <= 0);
  assert(abs(out.count(2) - 1) <= 0);

  fprintf('PASS: test_reconstruct_bddc_solution (reconstruct_bddc_solution)\n');
end

% ---------- local helpers ----------

function data = build_tiny_data_()
  % Two subdomains, each with one product interface dof scattered to the SAME global gamma dof.
  % Sub(1) also has one interior dof (scattered to global free dof #2).
  sub = struct([]);

  % --- subdomain 1 ---
  sub(1).prod_idx    = 1;
  sub(1).K_II        = 1;      % scalar (nI = 1)
  sub(1).R_II        = 1;      % chol factor of K_II (so K_II = R_II'*R_II)
  sub(1).f_I         = 10;     % scalar
  sub(1).K_Ig        = 1;      % (1x1) coupling to gamma value
  sub(1).glob_I      = 2;      % interior goes to global dof #2
  sub(1).gamma_glob  = 1;      % interface goes to global dof #1

  % --- subdomain 2 (no interior dofs) ---
  sub(2).prod_idx    = 2;
  sub(2).K_II        = zeros(0,0);
  sub(2).R_II        = 1;          % unused when nI == 0
  sub(2).f_I         = zeros(0,1);
  sub(2).K_Ig        = zeros(0,1);
  sub(2).glob_I      = zeros(0,1);
  sub(2).gamma_glob  = 1;          % interface also scatters to global dof #1 (averaging exercised)

  data = struct();
  data.sub = sub;

  % Hat-to-product trace operator: identity in the tiny test
  data.bddc = struct();
  data.bddc.R = speye(2);

  % Exactly two global free dofs: (1) shared gamma, (2) interior from sub(1)
  data.ddm = struct('nFree', 2);
end

function assert_throws_(fh, label)
  didThrow = false;
  try
    fh();
  catch
    didThrow = true;
  end
  if ~didThrow
    fprintf(2, 'Expected an error but none was thrown (%s)\n', label);
  end
  assert(didThrow);
end

function cleanup_paths_and_cache_(orig_path)
  % Restore path without prompting; clear function caches so persistent guards reset.
  try
    path(orig_path);
  catch
    % Best-effort; do not mask test failures with cleanup errors.
  end

  % Clear loaded functions (includes persistent state inside local helpers).
  try
    clear functions;
  catch
  end
end

function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did
    return;
  end

  % Strategy 1: setup_paths already resolvable
  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    main_dir = fileparts(sp);
    root_dir = fileparts(main_dir);
    addpath(genpath(fullfile(root_dir, 'tests')));
    did = true;
    return;
  end

  % Strategy 2: derive root from this test file location
  thisdir = fileparts(mfilename('fullpath'));
  root_dir = thisdir;
  found = false;

  for k = 1:10
    cand = fullfile(root_dir, 'main', 'setup_paths.m');
    if exist(cand, 'file') == 2
      found = true;
      break;
    end
    parent = fileparts(root_dir);
    if isempty(parent) || strcmp(parent, root_dir)
      break;
    end
    root_dir = parent;
  end

  if ~found
    error('could not locate project root containing main/setup_paths.m');
  end

  addpath(fullfile(root_dir, 'main'));
  setup_paths();
  addpath(genpath(fullfile(root_dir, 'tests')));

  did = true;
end