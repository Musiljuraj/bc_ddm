function test_apply_blockdiag_S()
%TEST_APPLY_BLOCKDIAG_S  Validate block-diagonal Schur apply in product space.
%
% Checks:
%   (1) apply_blockdiag_S(sub, w) matches explicit block-diagonal action
%       using the assembled local Schur complements S^(i).
%   (2) Basic input/consistency error paths (row vector, missing prod_idx,
%       out-of-range prod_idx) throw errors.

  fprintf('Running test_apply_blockdiag_S...\n');

  ensure_project_paths_();  % % CHANGED: use same path setup pattern as test_apply_local_schur.m

  tol = 1e-10;

  %------------------------------------------------------------
  % Build a small, deterministic pipeline
  %------------------------------------------------------------
  n = 6;
  nSubX = 3;
  nSubY = 2;

  [p, t, bnd] = mesh_unit_square_P1(n);
  f = @(x,y) 1 + x + y;

  [sub, ddm]   = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
  [sub, iface] = identify_interface_dofs(sub, ddm);

  sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f);
  sub = extract_subdomain_blocks(sub);

  opts = struct();
  opts.assemble_S = true;   % for explicit reference apply
  sub = setup_local_schur(sub, opts);

  [sub, prod] = build_product_interface(sub, iface);

  %------------------------------------------------------------
  % (1) Compare operator apply against explicit block diagonal.
  %------------------------------------------------------------
  w = (1:prod.nProd).';
  y1 = apply_blockdiag_S(sub, w);

  y2 = zeros(prod.nProd, 1);
  for i = 1:ddm.Nsub
    idx = sub(i).prod_idx(:);
    if isempty(idx)
      continue;
    end
    y2(idx) = sub(i).S * w(idx);
  end

  assert(numel(y1) == prod.nProd && size(y1,2) == 1, 'Output y has wrong shape.');
  assert(norm(full(y1) - full(y2), 2) < tol, ...
         'apply_blockdiag_S does not match explicit block-diagonal Schur action.');
  fprintf('  Explicit block-diagonal comparison: PASSED\n');

  %------------------------------------------------------------
  % (2) Error paths
  %------------------------------------------------------------

  % (2a) w must be a column vector
  threw = false;
  try
    apply_blockdiag_S(sub, w.');
  catch
    threw = true;
  end
  assert(threw, 'Expected error for row-vector input w.');
  fprintf('  Row-vector input check: PASSED\n');

  % (2b) missing prod_idx must error
  % IMPORTANT (Octave): cannot remove a field from only one element of a struct array.
  % Remove from the whole array; apply_blockdiag_S should fail immediately at i=1.
  threw = false;
  try
    sub_bad = rmfield(sub, 'prod_idx');
    apply_blockdiag_S(sub_bad, w);
  catch
    threw = true;
  end
  assert(threw, 'Expected error for missing sub(i).prod_idx.');
  fprintf('  Missing prod_idx check: PASSED\n');

  % (2c) out-of-range prod_idx must error
  sub_bad = sub;
  did_corrupt = false;
  for k = 1:numel(sub_bad)
    if ~isempty(sub_bad(k).prod_idx)
      sub_bad(k).prod_idx(1) = prod.nProd + 1;
      did_corrupt = true;
      break;
    end
  end
  assert(did_corrupt, 'Test setup error: could not find a non-empty prod_idx to corrupt.');

  threw = false;
  try
    apply_blockdiag_S(sub_bad, w);
  catch
    threw = true;
  end
  assert(threw, 'Expected error for out-of-range prod_idx.');
  fprintf('  Out-of-range prod_idx check: PASSED\n');

  fprintf('PASS: test_apply_blockdiag_S\n');
end

% =========================
% Helpers
% =========================

function ensure_project_paths_()  % % ADDED
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

% % FIXED: removed auto-execution "test_apply_blockdiag_S();" at EOF