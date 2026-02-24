function test_build_jump_operator_B()
% test_build_jump_operator_B  Unit tests for build_jump_operator_B
%
% Tests: build_jump_operator_B
% Focus: jump constraints structure (+/-1), row_hat consistency, rejection
%
% Suggested location: tests/ch3_ddm/test_build_jump_operator_B.m

  ensure_project_paths_();

  TARGET = 'build_jump_operator_B'; % ADDED
  WHAT   = 'jump constraints structure (+/-1), row_hat consistency, rejection'; % ADDED

  % -------------------------
  % 1) Wrong-arg-count checks
  % -------------------------
  assert_throws_(@() build_jump_operator_B(),           'wrong-arg-count: nargin==0'); % CHANGED
  assert_throws_(@() build_jump_operator_B(struct(),1), 'wrong-arg-count: nargin==2'); % CHANGED

  % -----------------------------------------
  % 2) Basic input validation (must error)
  % -----------------------------------------
  assert_throws_(@() build_jump_operator_B(123), 'invalid-type: prod not a struct'); % CHANGED

  s = struct();
  assert_throws_(@() build_jump_operator_B(s), 'missing-fields: hat2prod/nProd/nHat'); % CHANGED

  s.hat2prod = {1};
  s.nProd = 1;
  assert_throws_(@() build_jump_operator_B(s), 'missing-field: nHat'); % CHANGED

  % hat2prod present but not a cell -> should error (either explicit or from ''{ }'' indexing)
  s2 = struct();
  s2.hat2prod = [1 2 3];
  s2.nProd = 3;
  s2.nHat = 3;
  assert_throws_(@() build_jump_operator_B(s2), 'invalid-type: hat2prod not a cell'); % CHANGED

  % ---------------------------------------------------------
  % 3) Happy-path: small deterministic construction/invariants
  % ---------------------------------------------------------
  prod = struct();
  prod.nHat = 3;
  prod.nProd = 6;
  prod.hat2prod = { [1], [2 4], [3 5 6] };

  [B, meta] = build_jump_operator_B(prod);

  % Expected constraint count: (2-1) + (3-1) = 3
  assert(isequal(size(B), [3, 6]), 'unexpected size(B)'); % CHANGED
  assert(isstruct(meta) && isfield(meta, 'row_hat'), 'meta.row_hat missing'); % CHANGED
  assert(isequal(size(meta.row_hat), [3, 1]), 'unexpected size(meta.row_hat)'); % CHANGED

  % Entries must be in {-1, +1} for nonzeros
  nz = nonzeros(B);
  assert(~isempty(nz), 'B has no nonzeros but should'); % CHANGED
  assert(all(abs(nz) == 1), 'B nonzeros must be +/-1'); % CHANGED

  % Each constraint row should have exactly 2 nonzeros
  row_nnz = full(sum(B ~= 0, 2));
  assert(all(row_nnz == 2), 'each row must have exactly 2 nonzeros'); % CHANGED

  % Each row sums to zero ( +1 and -1 )
  row_sum = full(sum(B, 2));
  assert(all(row_sum == 0), 'each row must sum to zero'); % CHANGED

  % B * ones == 0
  z = full(B * ones(prod.nProd, 1));
  assert(norm(z, 2) == 0, 'B*ones must be exactly zero'); % CHANGED

  % Check per-row structure using meta.row_hat without assuming row ordering
  for r = 1:size(B,1)
    h = meta.row_hat(r);
    assert(h >= 1 && h <= prod.nHat, 'meta.row_hat out of range'); % CHANGED

    copies = prod.hat2prod{h}(:);
    assert(numel(copies) >= 2, 'row_hat must reference only hats with >=2 copies'); % CHANGED

    cols = find(B(r, :) ~= 0);
    vals = full(B(r, cols));

    % identify +1 and -1 columns
    p_plus  = cols(vals ==  1);
    p_minus = cols(vals == -1);

    assert(numel(p_plus) == 1,  'row must have exactly one +1'); % CHANGED
    assert(numel(p_minus) == 1, 'row must have exactly one -1'); % CHANGED

    p1 = copies(1);
    others = copies(2:end);

    assert(p_plus == p1, 'the +1 column must be copies(1)'); % CHANGED
    assert(any(p_minus == others), 'the -1 column must be one of copies(2:end)'); % CHANGED
  end

  % Check that for each hat, the set of "minus" columns equals copies(2:end)
  hats = unique(meta.row_hat(:)).';
  for h = hats
    copies = prod.hat2prod{h}(:);
    p1 = copies(1);
    others = copies(2:end);

    rows_h = find(meta.row_hat == h);
    minus_cols = zeros(numel(rows_h), 1);
    plus_cols  = zeros(numel(rows_h), 1);

    for ii = 1:numel(rows_h)
      r = rows_h(ii);
      cols = find(B(r, :) ~= 0);
      vals = full(B(r, cols));
      plus_cols(ii)  = cols(vals ==  1);
      minus_cols(ii) = cols(vals == -1);
    end

    % all plus columns must equal the reference p1
    assert(all(plus_cols == p1), 'all +1 columns must equal copies(1) for a hat'); % CHANGED

    % minus columns should match "others" as a set (ordering-independent)
    assert(isequal(sort(minus_cols), sort(others)), 'minus columns must equal copies(2:end) as a set'); % CHANGED
  end

  % -------------------------------------------
  % 4) Edge case: no constraints (m = 0) works
  % -------------------------------------------
  prod0 = struct();
  prod0.nHat = 2;
  prod0.nProd = 2;
  prod0.hat2prod = { [1], [2] };

  [B0, meta0] = build_jump_operator_B(prod0);
  assert(isequal(size(B0), [0, 2]), 'unexpected size(B0)'); % CHANGED
  assert(isstruct(meta0) && isfield(meta0, 'row_hat'), 'meta0.row_hat missing'); % CHANGED
  assert(isequal(size(meta0.row_hat), [0, 1]), 'unexpected size(meta0.row_hat)'); % CHANGED

  % -------------------------------------------------------
  % 5) Index-range / index-type: should error (from sparse)
  % -------------------------------------------------------
  prod_bad = struct();
  prod_bad.nHat = 1;
  prod_bad.nProd = 2;
  prod_bad.hat2prod = { [1 3] }; % 3 out of range for nProd=2
  assert_throws_(@() build_jump_operator_B(prod_bad), 'bad-index: out of range'); % CHANGED

  prod_bad2 = struct();
  prod_bad2.nHat = 1;
  prod_bad2.nProd = 3;
  prod_bad2.hat2prod = { [1 2.5] }; % non-integer index
  assert_throws_(@() build_jump_operator_B(prod_bad2), 'bad-index: non-integer'); % CHANGED

  prod_bad3 = struct();
  prod_bad3.nHat = 1;
  prod_bad3.nProd = 3;
  prod_bad3.hat2prod = { [1 NaN] }; % invalid index
  assert_throws_(@() build_jump_operator_B(prod_bad3), 'bad-index: NaN'); % CHANGED

  % ADDED: explicit PASS information when running this test directly
  fprintf('[PASS] %s  | tests: %s  | %s\n', mfilename(), TARGET, WHAT);
end

% -------------------------------------------------------------------------
function ensure_project_paths_()
% ensure_project_paths_  Idempotently add bc_ddm paths so tests run from anywhere.

  persistent did;
  if ~isempty(did) && did
    return;
  end

  % Strategy (1): setup_paths already resolvable
  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    main_dir = fileparts(sp);
    root_dir = fileparts(main_dir);
    addpath(genpath(fullfile(root_dir, 'tests')));
    did = true;
    return;
  end

  % Strategy (2): derive root from this test file location, walk upwards
  thisdir = fileparts(mfilename('fullpath'));
  root_dir = thisdir;

  found = false;
  for k = 1:10
    candidate = fullfile(root_dir, 'main', 'setup_paths.m');
    if exist(candidate, 'file') == 2
      found = true;
      break;
    end
    parent = fileparts(root_dir);
    if isempty(parent) || isequal(parent, root_dir)
      break;
    end
    root_dir = parent;
  end

  if ~found
    error('ensure_project_paths_: could not locate project root containing main/setup_paths.m');
  end

  addpath(fullfile(root_dir, 'main'));
  setup_paths();
  addpath(genpath(fullfile(root_dir, 'tests')));

  did = true;
end

% -------------------------------------------------------------------------
function assert_throws_(fh, label)
% assert_throws_  Assert that calling fh() throws an error (message not checked).
% label is used to make "expected error but none thrown" actionable.

  threw = false;
  try
    fh();
  catch
    threw = true;
  end

  if nargin < 2
    label = 'unlabeled-case';
  end

  assert(threw, sprintf('Expected an error but none was thrown (%s).', label)); % CHANGED
end