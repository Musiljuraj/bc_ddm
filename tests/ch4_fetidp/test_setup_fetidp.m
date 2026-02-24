function test_setup_fetidp()
%TEST_SETUP_FETIDP  Unit tests for setup_fetidp (FETI-DP setup).
%
% Suggested location: tests/ch4_feti_dp/test_setup_fetidp.m
%
% Covers:
%  - delta_idx_prod / delta_range mapping invariants
%  - local Schur split blocks Scc/Scd/Sdc/Sdd + Sdd Cholesky
%  - Bd construction + row reduction consistency vs Bd_full row subset
%  - DeltaWeights consistency vs multiplicity_scaling(prod)
%  - coarse operator Kcc/hc sizes + Cholesky (when nC>0)
%  - signature and missing-field error paths
%
% Path robustness:
%  - Uses ensure_project_paths_() matching the verified pattern in test_apply_local_schur.m

  ensure_project_paths_();
  seed_rng_();

  % -------------------------
  % Positive (success) path
  % -------------------------
  data = make_tiny_problem_data_();

  data_f = setup_fetidp(data, struct('tol_rank', 1e-12));

  % Required output fields
  req = {'delta_idx_prod','nDeltaProd','delta_range', ...
         'Scc','Scd','Sdc','Sdd','Sdd_R', ...
         'Bd','BdT','lambda_rows','nLambda', ...
         'DeltaWeights','Kcc','hc','Kcc_R'};
  for k = 1:numel(req)
    assert(isfield(data_f, req{k}), sprintf('Missing output field: %s', req{k}));
  end
  assert(isfield(data_f,'primal') && isfield(data_f.primal,'delta_idx'), 'data.primal.delta_idx missing');
  assert(isfield(data_f.primal,'prod2delta'), 'data.primal.prod2delta missing');
  assert(isfield(data_f.primal,'nC'), 'data.primal.nC missing');

  % delta_idx_prod must equal primal.delta_idx(:)
  assert(isequal(data_f.delta_idx_prod(:), data_f.primal.delta_idx(:)), 'delta_idx_prod mismatch');
  assert(data_f.nDeltaProd == numel(data_f.delta_idx_prod), 'nDeltaProd mismatch');

  % delta_range must match prod2delta(prod_idx_d) per subdomain
  nSub = numel(data_f.sub);
  assert(iscell(data_f.delta_range) && numel(data_f.delta_range) == nSub, 'delta_range shape mismatch');
  for i = 1:nSub
    assert(isfield(data_f.sub(i),'prod_idx_d'), 'sub(i).prod_idx_d missing');
    exp_range = data_f.primal.prod2delta(data_f.sub(i).prod_idx_d(:));
    got_range = data_f.delta_range{i}(:);
    assert(isequal(got_range, exp_range), sprintf('delta_range{%d} mapping mismatch', i));
  end

  % Local Schur split block sizes + Sdd factor consistency
  assert(numel(data_f.Scc)==nSub && numel(data_f.Sdd)==nSub, 'Schur cell sizes mismatch');
  for i = 1:nSub
    idx_c = data_f.sub(i).idx_c(:);
    idx_d = data_f.sub(i).idx_d(:);

    assert(all(size(data_f.Scc{i}) == [numel(idx_c), numel(idx_c)]), 'Scc size mismatch');
    assert(all(size(data_f.Scd{i}) == [numel(idx_c), numel(idx_d)]), 'Scd size mismatch');
    assert(all(size(data_f.Sdc{i}) == [numel(idx_d), numel(idx_c)]), 'Sdc size mismatch');
    assert(all(size(data_f.Sdd{i}) == [numel(idx_d), numel(idx_d)]), 'Sdd size mismatch');

    if isempty(idx_d)
      assert(isempty(data_f.Sdd_R{i}), 'Expected empty Sdd_R for empty idx_d');
    else
      R = data_f.Sdd_R{i};
      A = full(data_f.Sdd{i});
      assert_upper_triangular_(R, 1e-12);
      assert_close_(R' * R, A, 1e-10, 0, sprintf('Sdd factorization mismatch on subdomain %d', i));
    end
  end

  % Jump operator restriction: Bd must equal Bd_full(lambda_rows,:)
  assert(exist('build_jump_operator_B','file') == 2, 'build_jump_operator_B not found on path');
  B = build_jump_operator_B(data_f.prod);
  Bd_full = B(:, data_f.delta_idx_prod);

  assert(size(data_f.Bd,2) == size(Bd_full,2), 'Bd column count mismatch');
  assert_close_(full(data_f.BdT), full(data_f.Bd'), 0, 0, 'BdT must equal Bd'' exactly.');

  if ~isempty(Bd_full)
    lr = data_f.lambda_rows(:);
    assert(all(lr >= 1 & lr <= size(Bd_full,1)), 'lambda_rows out of range');
    assert_close_(full(data_f.Bd), full(Bd_full(lr,:)), 1e-12, 0, 'Bd must equal Bd_full(lambda_rows,:)');
    assert_full_row_rank_(data_f.Bd, 1e-12);
  end
  assert(data_f.nLambda == size(data_f.Bd,1), 'nLambda mismatch');

  % DeltaWeights must match multiplicity_scaling(prod) restricted to delta indices
  assert(exist('multiplicity_scaling','file') == 2, 'multiplicity_scaling not found on path');
  omega = multiplicity_scaling(data_f.prod);
  assert_close_(data_f.DeltaWeights(:), omega(data_f.delta_idx_prod(:)), 0, 0, 'DeltaWeights mismatch vs multiplicity_scaling(prod).');
  assert(all(isfinite(data_f.DeltaWeights(:))), 'DeltaWeights must be finite');
  assert(all(data_f.DeltaWeights(:) > 0), 'DeltaWeights must be positive');

  % Coarse operator checks
  nC = data_f.primal.nC;
  assert(all(size(data_f.Kcc) == [nC,nC]), 'Kcc size mismatch');
  assert(all(size(data_f.hc) == [nC,1]), 'hc size mismatch');

  if nC > 0
    assert_symmetric_(data_f.Kcc, 1e-12);
    assert_upper_triangular_(data_f.Kcc_R, 1e-12);
    assert_close_(data_f.Kcc_R' * data_f.Kcc_R, full(data_f.Kcc), 1e-10, 0, 'Kcc factorization mismatch.');
  else
    assert(isempty(data_f.Kcc_R), 'Expected empty Kcc_R when nC==0');
  end

  % -------------------------
  % Negative tests
  % -------------------------

  % Wrong arg count must throw
  assert_throws_(@() setup_fetidp(), 'setup_fetidp nargin=0');
  assert_throws_(@() setup_fetidp(struct(), struct(), 7), 'setup_fetidp nargin=3');

  % Missing sub(i).prod_idx_d must throw (delta_range construction)
  data_bad = struct();
  data_bad.sub = struct(); % 1 subdomain
  data_bad.primal = struct('delta_idx', 1, 'prod2delta', 1, 'nC', 0);
  data_bad.prod = struct();
  assert_throws_(@() setup_fetidp(data_bad), 'missing sub(i).prod_idx_d must error');

  % Missing sub(i).S must throw (local Schur requirement)
  data_bad2 = struct();
  data_bad2.sub = struct('prod_idx_d', 1); % no S
  data_bad2.primal = struct('delta_idx', 1, 'prod2delta', 1, 'nC', 0);
  data_bad2.prod = struct();
  assert_throws_(@() setup_fetidp(data_bad2), 'missing sub(i).S must error');

  % Non-SPD Sdd must throw (chol failure), before jump/scaling usage
  data_bad3 = struct();
  sub3 = struct();
  sub3.prod_idx_d = 1;
  sub3.idx_c = [];
  sub3.idx_d = 1;
  sub3.S = sparse(-1); % Sdd = [-1]
  sub3.g = 0;
  sub3.c_ids = [];
  data_bad3.sub = sub3;
  data_bad3.primal = struct('delta_idx', 1, 'prod2delta', 1, 'nC', 0);
  data_bad3.prod = struct();
  assert_throws_(@() setup_fetidp(data_bad3), 'non-SPD Sdd must error');

  fprintf('PASS: test_setup_fetidp\n');
end

% =========================
% Helpers
% =========================

function ensure_project_paths_()
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

function seed_rng_()
  if exist('rng','file') == 2
    rng(0);
  else
    rand('state', 0); %#ok<RAND>
    randn('state', 0); %#ok<RANDN>
  end
end

function data = make_tiny_problem_data_()
% Construct a small, deterministic "data" struct consumable by setup_fetidp.
%
% This avoids calling build_problem_data() with too few inputs (your repo
% uses narginchk), and instead tries a small set of common *positional*
% signatures. If build_problem_data returns a struct not yet "FETI-DP ready",
% we attempt to complete the pipeline inside the test (no project edits).

  if exist('build_problem_data','file') ~= 2
    error('test_setup_fetidp: build_problem_data not found on path.');
  end

  % Small-but-nontrivial parameters (interface + coarse should exist in typical setups)
  nSubX = 2;
  nSubY = 1;
  n = 12;
  f = @(x,y) 1 + 0*x;

  have_mesh = (exist('mesh_unit_square_P1','file') == 2);
  if have_mesh
    [p, t, bnd] = mesh_unit_square_P1(n);
  else
    p = []; t = []; bnd = [];
  end

  tried = {};
  lastErr = [];

  % Attempt 1: signature-driven call if we can parse it
  try
    [ok, data, label] = try_signature_driven_build_(p, t, bnd, n, nSubX, nSubY, f);
    if ok
      tried{end+1} = label;
      data = complete_pipeline_for_fetidp_(data);
      if is_valid_for_setup_fetidp_(data), return; end
    end
  catch err
    lastErr = err;
  end

  % Attempt 2: common positional patterns (mesh-based)
  if have_mesh
    try
      tried{end+1} = 'build_problem_data(p,t,bnd,nSubX,nSubY,f)';
      data = build_problem_data(p, t, bnd, nSubX, nSubY, f);
      data = complete_pipeline_for_fetidp_(data);
      if is_valid_for_setup_fetidp_(data), return; end
    catch err
      lastErr = err;
    end

    try
      tried{end+1} = 'build_problem_data(p,t,bnd,nSubX,nSubY,f,opts)';
      data = build_problem_data(p, t, bnd, nSubX, nSubY, f, struct());
      data = complete_pipeline_for_fetidp_(data);
      if is_valid_for_setup_fetidp_(data), return; end
    catch err
      lastErr = err;
    end
  end

  % Attempt 3: common positional patterns (grid-based)
  try
    tried{end+1} = 'build_problem_data(n,nSubX,nSubY,f)';
    data = build_problem_data(n, nSubX, nSubY, f);
    data = complete_pipeline_for_fetidp_(data);
    if is_valid_for_setup_fetidp_(data), return; end
  catch err
    lastErr = err;
  end

  try
    tried{end+1} = 'build_problem_data(n,nSubX,nSubY,f,opts)';
    data = build_problem_data(n, nSubX, nSubY, f, struct());
    data = complete_pipeline_for_fetidp_(data);
    if is_valid_for_setup_fetidp_(data), return; end
  catch err
    lastErr = err;
  end

  % Diagnostics: show what we tried and the likely fix
  fprintf('DIAG: test_setup_fetidp could not construct valid data for setup_fetidp.\n');
  fprintf('DIAG: Tried build_problem_data call patterns:\n');
  for k = 1:numel(tried)
    fprintf('  - %s\n', tried{k});
  end
  sp = which('build_problem_data');
  if ~isempty(sp)
    fprintf('DIAG: build_problem_data resolved to: %s\n', sp);
    sigline = parse_function_signature_line_(sp, 'build_problem_data');
    if ~isempty(sigline)
      fprintf('DIAG: Parsed signature line: %s\n', sigline);
    end
  end
  fprintf('DIAG: Next step: open your verified test_build_problem_data.m and mirror its call signature here.\n');

  if ~isempty(lastErr)
    rethrow(lastErr);
  end
  error('test_setup_fetidp: could not build a valid data struct for setup_fetidp.');
end

function [ok, data, label] = try_signature_driven_build_(p, t, bnd, n, nSubX, nSubY, f)
% Try to parse build_problem_data(...) signature and call it by argument names.
  ok = false; data = []; label = '';

  sp = which('build_problem_data');
  if isempty(sp)
    return;
  end
  args = parse_arglist_from_signature_(sp, 'build_problem_data');
  if isempty(args)
    return;
  end

  vals = cell(numel(args),1);
  for k = 1:numel(args)
    a = lower(strtrim(args{k}));

    if any(strcmp(a, {'p','points'}))
      vals{k} = p;
    elseif any(strcmp(a, {'t','tri','tris','elements'}))
      vals{k} = t;
    elseif any(strcmp(a, {'bnd','boundary','b'}))
      vals{k} = bnd;
    elseif any(strcmp(a, {'n','nref','nmesh'}))
      vals{k} = n;
    elseif any(strcmp(a, {'nsubx','nsx'}))
      vals{k} = nSubX;
    elseif any(strcmp(a, {'nsuby','nsy'}))
      vals{k} = nSubY;
    elseif strcmp(a, 'f') || strcmp(a, 'rhs') || strcmp(a, 'force')
      vals{k} = f;
    elseif strcmp(a, 'opts') || strcmp(a, 'options')
      vals{k} = struct();
    else
      % Unknown arg name => abort signature-driven attempt (avoid guessing wrong types)
      return;
    end
  end

  label = sprintf('signature-driven build_problem_data(%s)', strjoin(args, ','));
  data = build_problem_data(vals{:});
  ok = true;
end

function data = complete_pipeline_for_fetidp_(data)
% If build_problem_data returns a struct missing pieces setup_fetidp expects,
% try to complete those pieces via known pipeline stages (best-effort).

  if ~isstruct(data) || ~isfield(data,'sub')
    return;
  end

  % Ensure explicit local Schur S and condensed RHS g exist (assemble_S=true)
  try
    if exist('setup_local_schur','file') == 2
      needsS = false;
      for i = 1:numel(data.sub)
        if ~isfield(data.sub(i),'S') || isempty(data.sub(i).S)
          needsS = true; break;
        end
      end
      if needsS
        data.sub = setup_local_schur(data.sub, struct('assemble_S', true));
      end
    end
  catch
    % best-effort only; leave as-is
  end

  % Ensure primal maps: prod_idx_d, idx_c, idx_d, c_ids (best-effort)
  try
    if exist('build_primal_maps','file') == 2
      needsPrimalMaps = false;
      for i = 1:numel(data.sub)
        si = data.sub(i);
        if ~isfield(si,'prod_idx_d') || ~isfield(si,'idx_c') || ~isfield(si,'idx_d') || ~isfield(si,'c_ids')
          needsPrimalMaps = true; break;
        end
      end
      if needsPrimalMaps
        % Common pattern is build_primal_maps(data) -> data
        data = build_primal_maps(data);
      end
    end
  catch
    % best-effort only; leave as-is
  end
end

function tf = is_valid_for_setup_fetidp_(data)
% Strong “ready for setup_fetidp” validator.

  tf = isstruct(data) && isfield(data,'sub') && isfield(data,'primal') && isfield(data,'prod');
  if ~tf, return; end

  tf = tf ...
    && isfield(data.primal,'delta_idx') ...
    && isfield(data.primal,'prod2delta') ...
    && isfield(data.primal,'nC');
  if ~tf, return; end

  if ~isstruct(data.sub) || isempty(data.sub)
    tf = false; return;
  end

  % At least one subdomain must contain all key fields used by setup_fetidp.
  okOne = false;
  for i = 1:numel(data.sub)
    si = data.sub(i);
    ok = isfield(si,'prod_idx_d') && isfield(si,'idx_c') && isfield(si,'idx_d') && isfield(si,'c_ids') ...
         && isfield(si,'S') && ~isempty(si.S) ...
         && isfield(si,'g');
    if ok
      okOne = true;
      break;
    end
  end
  tf = okOne;
end

function sigline = parse_function_signature_line_(filepath, fname)
% Return the first "function ..." line mentioning fname, else ''.
  sigline = '';
  try
    txt = fileread(filepath);
  catch
    return;
  end
  lines = regexp(txt, '\r\n|\n|\r', 'split');
  for i = 1:min(60, numel(lines))
    L = strtrim(lines{i});
    if startsWith(L, 'function') && ~isempty(strfind(L, fname))
      sigline = L;
      return;
    end
  end
end

function args = parse_arglist_from_signature_(filepath, fname)
% Parse "function ... = fname(arg1,arg2,...)" and return {arg1,arg2,...}.
  args = {};
  sigline = parse_function_signature_line_(filepath, fname);
  if isempty(sigline)
    return;
  end
  % Match fname(...)
  pat = [fname '\s*\(([^)]*)\)'];
  tok = regexp(sigline, pat, 'tokens', 'once');
  if isempty(tok)
    return;
  end
  inside = strtrim(tok{1});
  if isempty(inside)
    args = {};
    return;
  end
  parts = regexp(inside, ',', 'split');
  args = cellfun(@(s) strtrim(s), parts, 'UniformOutput', false);
end

function assert_upper_triangular_(R, tol)
  if nargin < 2, tol = 0; end
  L = tril(R, -1);
  ok = norm(L(:), 2) <= tol * max(1, norm(R(:),2));
  assert(ok, 'Expected (numerically) upper triangular factor.');
end

function assert_symmetric_(A, tol)
  if nargin < 2, tol = 0; end
  Af = full(A);
  assert_close_(Af, Af', tol, 0, 'Expected symmetric matrix.');
end

function assert_full_row_rank_(A, reltol)
  if nargin < 2, reltol = 1e-12; end
  Af = full(A);
  if isempty(Af)
    return;
  end
  s = svd(Af);
  if isempty(s)
    return;
  end
  r = sum(s > reltol * max(s));
  assert(r == size(Af,1), 'Expected full row rank (numerical).');
end

function assert_close_(a, b, rtol, atol, msg)
  if nargin < 3, rtol = 0; end
  if nargin < 4, atol = 0; end
  if nargin < 5, msg = 'Values not close.'; end

  da = full(a); db = full(b);
  err = norm(da - db, inf);
  denom = max(1, norm(db, inf));
  ok = (err <= atol + rtol * denom);

  assert(ok, sprintf('%s (err=%g, denom=%g, rtol=%g, atol=%g)', msg, err, denom, rtol, atol));
end

function assert_throws_(fh, label)
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