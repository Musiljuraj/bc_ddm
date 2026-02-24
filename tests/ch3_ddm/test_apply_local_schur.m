function test_apply_local_schur()
%TEST_APPLY_LOCAL_SCHUR  Unit tests for apply_local_schur (matrix-free local Schur action).
%
% Covers:
%   - nI>0 branch (uses R_II backsolves) vs explicit Schur action
%   - nI==0 branch (reduces to K_gg*x)
%   - signature/shape checks + missing-field checks
%   - CONTRACT check: output of extract_subdomain_blocks + setup_local_schur is consumable

  ensure_project_paths_();

  % ----------------------------
  % Case A: nI > 0 (R_II nonempty) on synthetic fixture
  % ----------------------------
  [si, K_II] = make_fixture_si_interior_();
  nG = size(si.K_gg, 1);

  rng(0);
  x = randn(nG,1);
  y = apply_local_schur(si, x);

  % Reference via explicit Schur complement using K_II (kept only in the test)
  S_ref = si.K_gg - si.K_gI * (K_II \ si.K_Ig);
  y_ref = S_ref * x;

  assert_close_(y, y_ref, 500*eps, 0, 'apply_local_schur mismatch vs explicit Schur action (synthetic, nI>0).');

  % ----------------------------
  % Case B: nI == 0 (R_II empty => y = K_gg*x) on synthetic fixture
  % ----------------------------
  sj = make_fixture_si_no_interior_();
  x2 = [3; -1];
  y2 = apply_local_schur(sj, x2);
  assert_close_(y2, sj.K_gg * x2, 0, 0, 'apply_local_schur should reduce to K_gg*x when nI==0.');

  % ----------------------------
  % Case C: CONTRACT test with real pipeline (extract_subdomain_blocks + setup_local_schur)
  % ----------------------------
  si_real = make_real_pipeline_si_();     % returns a subdomain struct with fields + explicit S
  nG3 = size(si_real.K_gg, 1);

  rng(1);
  x3 = randn(nG3,1);
  y3 = apply_local_schur(si_real, x3);
  y3_ref = si_real.S * x3;

  assert_close_(y3, y3_ref, 1e-10, 0, 'apply_local_schur mismatch vs explicit Schur action (pipeline contract).');

  % ----------------------------
  % Signature / contract checks
  % ----------------------------
  assert_throws_(@() apply_local_schur(si), 'apply_local_schur nargin=1');
  assert_throws_(@() apply_local_schur(si, x, 7), 'apply_local_schur nargin=3');

  % x must be a column vector of length nG
  assert_throws_(@() apply_local_schur(si, x.'), 'row-vector x must error');
  assert_throws_(@() apply_local_schur(si, [1;2;3]), 'wrong-length x must error');

  % required fields
  si_bad = rmfield(si, 'R_II');
  assert_throws_(@() apply_local_schur(si_bad, x), 'missing R_II must error');

  si_bad = rmfield(si, 'K_gI');
  assert_throws_(@() apply_local_schur(si_bad, x), 'missing K_gI must error');

  si_bad = rmfield(si, 'K_Ig');
  assert_throws_(@() apply_local_schur(si_bad, x), 'missing K_Ig must error');

  si_bad = rmfield(si, 'K_gg');
  assert_throws_(@() apply_local_schur(si_bad, x), 'missing K_gg must error');

  fprintf('PASS: test_apply_local_schur\n');
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

function [si, K_II] = make_fixture_si_interior_()
% Minimal si struct with interior (R_II nonempty).
  K_II = [2, -1; -1, 2];   % SPD
  K_Ig = [1, 0; 2, 1];     % (2x2)
  K_gI = K_Ig.';           % symmetric-consistent
  K_gg = [3, 1; 1, 2];

  R_II = chol(K_II);
  si = struct('K_gg', K_gg, 'K_Ig', K_Ig, 'K_gI', K_gI, 'R_II', R_II);
end

function si = make_fixture_si_no_interior_()
% nI==0 case uses empty R_II; then operator reduces to K_gg.
  K_gg = [4, 1; 1, 2];
  K_Ig = zeros(0,2);
  K_gI = zeros(2,0);
  R_II = [];

  si = struct('K_gg', K_gg, 'K_Ig', K_Ig, 'K_gI', K_gI, 'R_II', R_II);
end

function si = make_real_pipeline_si_()
% Build one subdomain struct from the real DDM pipeline and ensure it has:
%   K_gg, K_Ig, K_gI, R_II, and explicit S (assemble_S=true).
%
% The parameters are chosen to reliably produce subdomains with nI>0 and nG>0.

  n = 12;
  nSubX = 2;
  nSubY = 2;
  f = @(x,y) 1 + x - 2*y;

  [p, t, bnd] = mesh_unit_square_P1(n);

  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
  [sub, iface] = identify_interface_dofs(sub, ddm); %#ok<NASGU>

  sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f);
  sub = extract_subdomain_blocks(sub);
  sub = setup_local_schur(sub, struct('assemble_S', true));

  % Pick a subdomain with nonempty interface AND nonempty interior so we exercise the R_II branch.
  picked = false;
  for i = 1:ddm.Nsub
    nG = size(sub(i).K_gg, 1);
    if nG > 0 && ~isempty(sub(i).R_II)
      si = sub(i);
      picked = true;
      break;
    end
  end

  % Fallback (should not trigger with the parameters above), but keep deterministic.
  if ~picked
    for i = 1:ddm.Nsub
      nG = size(sub(i).K_gg, 1);
      if nG > 0
        si = sub(i);
        return;
      end
    end
    error('test_apply_local_schur: could not find any subdomain with nG>0 for contract test.');
  end
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