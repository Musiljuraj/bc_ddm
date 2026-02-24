function test_reconstruct_fetidp_solution()
% Minimal-standard test for reconstruct_fetidp_solution()

  fprintf('Running test_reconstruct_fetidp_solution (function: reconstruct_fetidp_solution)\n');
  ensure_project_paths_();

  % --- wrong-arg-count must throw ---
  assert_throws_(@() reconstruct_fetidp_solution(), 'no-args');
  assert_throws_(@() reconstruct_fetidp_solution(1), 'one-arg');
  assert_throws_(@() reconstruct_fetidp_solution(1, struct(), 3), 'three-args');

  % --- representative invalid input must throw (dimension mismatch) ---
  data = build_tiny_data_();
  assert_throws_(@() reconstruct_fetidp_solution([1;2;3], data), 'lambda-length-mismatch');

  % --- representative invalid input must throw (missing required field) ---
  data2 = rmfield(data, 'BdT');
  assert_throws_(@() reconstruct_fetidp_solution([1;2], data2), 'missing-field-BdT');

  % --- smoke path + output contract checks ---
  lambda = [1; 2];
  [w_c, w_d, u_free, diag] = reconstruct_fetidp_solution(lambda, data);

  % Shapes / sizes
  assert(isvector(w_c));
  assert(numel(w_c) == data.primal.nC);

  assert(isvector(w_d));
  assert(numel(w_d) == data.nDeltaProd);

  assert(isvector(u_free));
  assert(numel(u_free) == numel(data.free));

  assert(isstruct(diag));
  assert(isfield(diag, 'constraint_norm'));

  % Real/finite numeric values where expected
  assert(isreal(w_d) && all(isfinite(w_d)));
  assert(isreal(u_free) && all(isfinite(u_free)));
  assert(isreal(diag.constraint_norm) && isfinite(diag.constraint_norm));

  % Stable invariant: diag.constraint_norm matches its defining expression
  cn = norm(full(data.Bd * w_d), 2);
  assert(abs(diag.constraint_norm - cn) <= 1e-12 * max(1, cn));

  % Stable behavior check: duplicate scatter averaging (both subs write to same global dof)
  % With this tiny construction, expected u_free is the average of the two local gamma values.
  assert(abs(u_free(1) - 2.5) <= 1e-12);

  fprintf('PASS: test_reconstruct_fetidp_solution (reconstruct_fetidp_solution)\n');
end

% ---------- local helpers ----------

function data = build_tiny_data_()
  % Two subdomains, each with one delta dof, both scattered onto the same global free dof.
  sub = struct([]);

  sub(1).g      = 3;
  sub(1).idx_c  = [];
  sub(1).idx_d  = 1;
  sub(1).c_ids  = [];
  sub(1).glob_G = 1;
  sub(1).glob_I = [];

  sub(2).g      = 5;
  sub(2).idx_c  = [];
  sub(2).idx_d  = 1;
  sub(2).c_ids  = [];
  sub(2).glob_G = 1;
  sub(2).glob_I = [];

  data = struct();
  data.sub = sub;
  data.nDeltaProd = 2;
  data.delta_range = {1, 2};

  % Identity BdT so t = lambda.
  data.BdT = speye(2);

  % Sdd_R factors are scalar 1 => Sdd^{-1} is identity in the tiny case.
  data.Sdd_R = {1, 1};

  % No coarse space in this tiny test.
  data.primal = struct('nC', 0);
  data.hc = zeros(0,1);
  data.Kcc_R = 1; % unused when nC=0, but set defensively.

  % Couplings not used because nC=0, but provide consistent empty shapes.
  data.Scd = {zeros(0,1), zeros(0,1)};
  data.Sdc = {zeros(1,0), zeros(1,0)};

  % One global FREE dof (both subs scatter to index 1 => averaging exercised).
  data.free = (1:1)';

  % Simple constraint operator (1 x 2).
  data.Bd = sparse(1, [1 2], [1 -1], 1, 2);
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