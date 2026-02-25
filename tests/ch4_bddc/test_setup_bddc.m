function test_setup_bddc()
%TEST_SETUP_BDDC  Unit tests for setup_bddc (BDDC setup).
%
% Suggested location: tests/ch4_bddc/test_setup_bddc.m
%
% Covers:
%  - signature / missing-field error paths
%  - g_prod assembly invariant
%  - omega == multiplicity_scaling(prod)
%  - R == build_assembly_operator_R(prod)
%  - local Sdd Cholesky factors Sdd_R
%  - coarse K0 / K0_R sizes + factorization (when nC>0)
%
% Path robustness:
%  - Uses ensure_project_paths_() matching the verified pattern in test_apply_local_schur.m

  ensure_project_paths_();
  seed_rng_();

  % -------------------------
  % Positive (success) path
  % -------------------------
  data = make_tiny_problem_data_();
  data = complete_pipeline_for_bddc_(data);

  assert(is_valid_for_setup_bddc_(data), 'Could not build valid data for setup_bddc.');

  data_b = setup_bddc(data, struct('tol_chol', 0));

  assert(isfield(data_b,'bddc') && isstruct(data_b.bddc), 'data.bddc missing or not a struct');
  b = data_b.bddc;

  % Required output fields
  req = {'R','omega','g_prod','Sdd_R','Psi','K0','K0_R'};
  for k = 1:numel(req)
    assert(isfield(b, req{k}), sprintf('Missing output field: data.bddc.%s', req{k}));
  end

  % Basic sizes
  assert(issparse(b.R), 'bddc.R must be sparse');
  assert(isvector(b.omega) && numel(b.omega) == data_b.prod.nProd, 'omega size mismatch');
  assert(isvector(b.g_prod) && numel(b.g_prod) == data_b.prod.nProd, 'g_prod size mismatch');
  assert(all(size(b.R) == [data_b.prod.nProd, data_b.prod.nHat]), 'R size mismatch');
  assert(all(size(b.Psi) == [data_b.prod.nProd, data_b.primal.nC]), 'Psi size mismatch');
  assert(all(size(b.K0)  == [data_b.primal.nC, data_b.primal.nC]), 'K0 size mismatch');

  % Finite checks (cheap)
  assert(all(isfinite(b.omega(:))), 'omega contains NaN/Inf');
  assert(all(isfinite(b.g_prod(:))), 'g_prod contains NaN/Inf');
  assert(all(isfinite(full(b.Psi(:)))), 'Psi contains NaN/Inf');
  assert(all(isfinite(b.K0(:))), 'K0 contains NaN/Inf');

  % Invariant: omega must match multiplicity_scaling(prod)
  assert(exist('multiplicity_scaling','file') == 2, 'multiplicity_scaling not found on path');
  omega_ref = multiplicity_scaling(data_b.prod);
  assert_close_(b.omega(:), omega_ref(:), 0, 0, 'omega mismatch vs multiplicity_scaling(prod).');
  assert(all(b.omega(:) > 0), 'omega must be positive');

  % Invariant: R must match build_assembly_operator_R(prod)
  assert(exist('build_assembly_operator_R','file') == 2, 'build_assembly_operator_R not found on path');
  R_ref = build_assembly_operator_R(data_b.prod);
  assert_close_(full(b.R), full(R_ref), 0, 0, 'R mismatch vs build_assembly_operator_R(prod).');

  % Invariant: g_prod assembled by scattering sub(i).g into prod_idx
  g_ref = zeros(data_b.prod.nProd, 1);
  for i = 1:numel(data_b.sub)
    idx = data_b.sub(i).prod_idx(:);
    if isempty(idx), continue; end
    g_ref(idx) = data_b.sub(i).g(:);
  end
  assert_close_(b.g_prod(:), g_ref(:), 0, 0, 'g_prod assembly mismatch.');

  % Local Sdd factorization checks
  nSub = numel(data_b.sub);
  assert(iscell(b.Sdd_R) && numel(b.Sdd_R) == nSub, 'Sdd_R must be cell(nSub,1)');
  for i = 1:nSub
    idx_d = data_b.sub(i).idx_d(:);
    if isempty(idx_d)
      assert(isempty(b.Sdd_R{i}), 'Expected empty Sdd_R for empty idx_d');
    else
      Rdd = b.Sdd_R{i};
      Sdd = full(data_b.sub(i).S(idx_d, idx_d));
      assert_upper_triangular_(Rdd, 1e-12);
      assert_close_(Rdd' * Rdd, Sdd, 1e-10, 0, sprintf('Sdd factorization mismatch on subdomain %d', i));
    end
  end

  % Coarse operator checks
  nC = data_b.primal.nC;
  if nC > 0
    assert_symmetric_(b.K0, 1e-12);
    assert_upper_triangular_(b.K0_R, 1e-12);
    assert_close_(b.K0_R' * b.K0_R, full(b.K0), 1e-10, 0, 'K0 factorization mismatch.');
  else
    assert(isempty(b.K0_R), 'Expected empty K0_R when nC==0');
    assert(isempty(b.K0) || all(size(b.K0)==[0,0]), 'Expected empty K0 when nC==0');
  end

  % -------------------------
  % Negative tests
  % -------------------------

  % Wrong arg count must throw
  assert_throws_(@() setup_bddc(), 'setup_bddc nargin=0');
  assert_throws_(@() setup_bddc(struct(), struct(), 7), 'setup_bddc nargin=3');

  % Missing top-level fields must throw
  assert_throws_(@() setup_bddc(struct('sub',[])), 'missing prod/primal must error');
  assert_throws_(@() setup_bddc(struct('prod',struct())), 'missing sub/primal must error');

  % Missing sub(i).g must throw
  data_bad = data;
  if isfield(data_bad.sub,'g')
    data_bad.sub = rmfield(data_bad.sub, 'g');
  end
  assert_throws_(@() setup_bddc(data_bad), 'missing sub(i).g must error');

  % Missing sub(i).S must throw
  data_bad2 = data;
  if isfield(data_bad2.sub,'S')
    data_bad2.sub = rmfield(data_bad2.sub, 'S');
  end
  assert_throws_(@() setup_bddc(data_bad2), 'missing sub(i).S must error');

  % Non-SPD Sdd must throw (chol failure)
  data_bad3 = data;
  [iSPD, idx_d] = find_sub_with_delta_(data_bad3.sub);
  if iSPD > 0
    S = data_bad3.sub(iSPD).S;
    m = numel(idx_d);
    S(idx_d, idx_d) = -speye(m); % force non-SPD
    data_bad3.sub(iSPD).S = S;
    assert_throws_(@() setup_bddc(data_bad3), 'non-SPD Sdd must error');
  end

  fprintf('PASS: test_setup_bddc\n');
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

function [iFound, idx_d] = find_sub_with_delta_(sub)
  iFound = 0; idx_d = [];
  for i = 1:numel(sub)
    if isfield(sub(i),'idx_d')
      id = sub(i).idx_d(:);
      if ~isempty(id)
        iFound = i;
        idx_d = id;
        return;
      end
    end
  end
end

function data = complete_pipeline_for_bddc_(data)
% Ensure setup_bddc prerequisites exist: sub(i).S, sub(i).g, idx_c/idx_d/c_ids.

  if ~isstruct(data) || ~isfield(data,'sub') || ~isfield(data,'prod') || ~isfield(data,'primal')
    return;
  end

  % Ensure local Schur provides explicit S and condensed g (assemble_S=true)
  try
    if exist('setup_local_schur','file') == 2
      needs = false;
      for i = 1:numel(data.sub)
        if ~isfield(data.sub(i),'S') || isempty(data.sub(i).S) || ~isfield(data.sub(i),'g')
          needs = true; break;
        end
      end
      if needs
        data.sub = setup_local_schur(data.sub, struct('assemble_S', true));
      end
    end
  catch
    % best-effort only
  end

  % Ensure primal maps exist (idx_c, idx_d, c_ids)
  try
    if exist('build_primal_maps','file') == 2
      needs = false;
      for i = 1:numel(data.sub)
        si = data.sub(i);
        if ~isfield(si,'idx_c') || ~isfield(si,'idx_d') || ~isfield(si,'c_ids')
          needs = true; break;
        end
      end
      if needs
        data = build_primal_maps(data);
      end
    end
  catch
    % best-effort only
  end
end

function tf = is_valid_for_setup_bddc_(data)
% Strong “ready for setup_bddc” validator.
  tf = isstruct(data) && isfield(data,'sub') && isfield(data,'prod') && isfield(data,'primal');
  if ~tf, return; end
  tf = tf && isfield(data.primal,'nC');
  if ~tf, return; end
  if ~isstruct(data.sub) || isempty(data.sub)
    tf = false; return;
  end

  okOne = false;
  for i = 1:numel(data.sub)
    si = data.sub(i);
    ok = isfield(si,'prod_idx') ...
         && isfield(si,'S') && ~isempty(si.S) ...
         && isfield(si,'g') ...
         && isfield(si,'idx_c') && isfield(si,'idx_d') && isfield(si,'c_ids');
    if ok
      okOne = true;
      break;
    end
  end
  tf = okOne;
end

function data = make_tiny_problem_data_()
% Construct a small, deterministic "data" struct consumable by setup_bddc.
% Mirrors the robust pattern used by test_setup_fetidp (tries several call signatures).

  if exist('build_problem_data','file') ~= 2
    error('test_setup_bddc: build_problem_data not found on path.');
  end

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
      if isstruct(data), return; end
    end
  catch err
    lastErr = err;
  end

  % Attempt 2: common positional patterns (mesh-based)
  if have_mesh
    try
      tried{end+1} = 'build_problem_data(p,t,bnd,nSubX,nSubY,f)';
      data = build_problem_data(p, t, bnd, nSubX, nSubY, f);
      if isstruct(data), return; end
    catch err
      lastErr = err;
    end

    try
      tried{end+1} = 'build_problem_data(p,t,bnd,nSubX,nSubY,f,opts)';
      data = build_problem_data(p, t, bnd, nSubX, nSubY, f, struct());
      if isstruct(data), return; end
    catch err
      lastErr = err;
    end
  end

  % Attempt 3: common positional patterns (grid-based)
  try
    tried{end+1} = 'build_problem_data(n,nSubX,nSubY,f)';
    data = build_problem_data(n, nSubX, nSubY, f);
    if isstruct(data), return; end
  catch err
    lastErr = err;
  end

  try
    tried{end+1} = 'build_problem_data(n,nSubX,nSubY,f,opts)';
    data = build_problem_data(n, nSubX, nSubY, f, struct());
    if isstruct(data), return; end
  catch err
    lastErr = err;
  end

  % Diagnostics
  fprintf('DIAG: test_setup_bddc could not construct data via build_problem_data.\n');
  fprintf('DIAG: Tried call patterns:\n');
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
  fprintf('DIAG: Next step: mirror the call signature used in your verified test_build_problem_data.m\n');

  if ~isempty(lastErr)
    rethrow(lastErr);
  end
  error('test_setup_bddc: could not build a data struct.');
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
      return; % unknown arg => abort (avoid guessing)
    end
  end

  label = sprintf('signature-driven build_problem_data(%s)', strjoin(args, ','));
  data = build_problem_data(vals{:});
  ok = true;
end

function sigline = parse_function_signature_line_(filepath, fname)
% Return the first "function ..." line mentioning fname, else ''.
  sigline = '';
  try
    txt = fileread(filepath);
  catch
    return;
  end
  lines = regexp(txt, '\r\n|\r|\n', 'split');
  for i = 1:numel(lines)
    s = strtrim(lines{i});
    if numel(s) >= 8 && strcmp(s(1:8), 'function')
      if ~isempty(strfind(s, fname)) %#ok<STREMP>
        sigline = s;
        return;
      end
    end
  end
end

function args = parse_arglist_from_signature_(filepath, fname)
% Parse argument list from the first function line for fname: function ... fname(a,b,c)
  args = {};
  sigline = parse_function_signature_line_(filepath, fname);
  if isempty(sigline)
    return;
  end
  pat = [fname '\s*\(([^)]*)\)'];
  tok = regexp(sigline, pat, 'tokens', 'once');
  if isempty(tok)
    return;
  end
  inside = strtrim(tok{1});
  if isempty(inside)
    return;
  end
  parts = regexp(inside, ',', 'split');
  args = cellfun(@(s) strtrim(s), parts, 'UniformOutput', false);
end

% -------------------------
% Assertions (mirrors verified helpers used in test_setup_fetidp)
% -------------------------

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
