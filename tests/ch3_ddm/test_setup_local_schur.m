function test_setup_local_schur()
%TEST_SETUP_LOCAL_SCHUR  Unit tests for setup_local_schur (DDM preparation stage).
%
% % CHANGED: split out from former integration-style combined test.
% % CHANGED: use tiny deterministic synthetic fixtures (no FEM/partition deps).
% % ADDED: path-robust helper (setup_paths + tests/ genpath).
% % ADDED: cover nI==0 branch, default opts behavior, multi-subdomain loop.
% % ADDED: signature checks + missing-field error checks.
% % ADDED: strict invalid-input rejection tests (char/logical/NaN/Inf/complex).
%   If any strict rejection test fails with "Expected an error but none was thrown",
%   treat as Rule 7.1 and harden setup_local_schur (do not relax the test).

  ensure_project_paths_();  % % ADDED

  % ----------------------------
  % Case A: multi-subdomain array; one with interior, one without
  % ----------------------------
  sub = repmat(struct(), 2, 1);  % % ADDED
  sub(1) = make_fixture_interior_();    % nI=2, nG=2
  sub(2) = make_fixture_no_interior_(); % nI=0, nG=2

  opts = struct('assemble_S', true);
  out  = setup_local_schur(sub, opts);

  % --- sub(1): nI>0
  s1 = out(1);
  assert(isfield(s1,'R_II') && isfield(s1,'g') && isfield(s1,'S'), ...
         'Expected R_II, g, S fields for nI>0 with assemble_S=true.');
  assert(isequal(size(s1.R_II), [2,2]), 'Expected R_II to be 2x2 for fixture nI=2.');
  assert(istriu(s1.R_II), 'Expected R_II to be upper triangular (chol factor).');
  assert(all(diag(s1.R_II) > 0), 'Expected positive diagonal in Cholesky factor R_II.');
  assert(isequal(size(s1.g), [2,1]), 'Expected g to be nGx1 (here 2x1).');
  assert(isequal(size(s1.S), [2,2]), 'Expected S to be nGxnG (here 2x2).');

  % g = f_g - K_gI*(K_II^{-1} f_I)
  g_ref = s1.f_g - s1.K_gI * (s1.K_II \ s1.f_I);
  assert_close_(s1.g, g_ref, 200*eps, 0, 'Reduced RHS g mismatch for nI>0.');

  % S = K_gg - K_gI*(K_II^{-1} K_Ig)
  S_ref = s1.K_gg - s1.K_gI * (s1.K_II \ s1.K_Ig);
  assert_close_(s1.S, S_ref, 500*eps, 0, 'Explicit Schur S mismatch for nI>0.');

  % For symmetric-consistent fixture, Schur should be symmetric (numerically)
  assert_close_(S_ref, S_ref.', 1e3*eps, 0, 'Schur complement should be numerically symmetric for symmetric blocks.');

  % --- sub(2): nI==0
  s2 = out(2);
  assert(isfield(s2,'R_II') && isfield(s2,'g') && isfield(s2,'S'), ...
         'Expected R_II, g, S fields for nI==0 with assemble_S=true.');
  assert(isempty(s2.R_II), 'Expected empty R_II when nI==0.');
  assert_close_(s2.g, s2.f_g, 0, 0, 'Expected g == f_g when nI==0.');
  assert_close_(s2.S, s2.K_gg, 0, 0, 'Expected S == K_gg when nI==0.');

  % ----------------------------
  % Case B: default opts (assemble_S defaults false => no S field)
  % ----------------------------
  out2 = setup_local_schur(sub); % nargin==1
  assert(isfield(out2(1),'R_II') && isfield(out2(1),'g'), 'Expected R_II and g fields with default opts.');
  assert(~isfield(out2(1),'S'), 'Did not expect field S when opts.assemble_S is default/false.');
  assert(~isfield(out2(2),'S'), 'Did not expect field S for nI==0 when opts.assemble_S is default/false.');

  % ----------------------------
  % Signature / error behavior
  % ----------------------------
  assert_throws_(@() setup_local_schur(), 'setup_local_schur nargin=0');                 % % ADDED
  assert_throws_(@() setup_local_schur(sub, opts, 7), 'setup_local_schur nargin=3');     % % ADDED

  % Octave struct arrays must have identical fields across all elements.
  % So for "missing field" checks, test on a scalar subdomain.
  sub_missing = sub(1);                                 % % FIXED
  sub_missing = rmfield(sub_missing, 'f_g');            % % FIXED
  assert_throws_(@() setup_local_schur(sub_missing, opts), 'missing f_g should error'); % % CHANGED

  sub_missing2 = sub(1);                                % % FIXED
  sub_missing2 = rmfield(sub_missing2, 'K_II');         % % FIXED
  assert_throws_(@() setup_local_schur(sub_missing2, opts), 'missing K_II should error'); % % CHANGED

  % ----------------------------
  % Strict invalid-input rejection (Rule 7.1)
  % ----------------------------
  % If any of these do NOT error, harden setup_local_schur input validation
  % (reject char/logical; require numeric real finite; enforce block sizes).

  sub_bad = make_fixture_interior_();
  sub_bad.K_II = char([2 1; 1 2]);  % char but SPD numerically -> must be rejected explicitly
  assert_throws_(@() setup_local_schur(sub_bad, opts), 'K_II is char');  % % ADDED

  sub_bad = make_fixture_interior_();
  sub_bad.K_II = logical([1 0; 0 1]); % logical SPD numerically -> must be rejected explicitly
  assert_throws_(@() setup_local_schur(sub_bad, opts), 'K_II is logical'); % % ADDED

  sub_bad = make_fixture_interior_();
  sub_bad.f_I = [NaN; 0];
  assert_throws_(@() setup_local_schur(sub_bad, opts), 'f_I contains NaN'); % % ADDED

  sub_bad = make_fixture_interior_();
  sub_bad.f_g = [Inf; 0];
  assert_throws_(@() setup_local_schur(sub_bad, opts), 'f_g contains Inf'); % % ADDED

  sub_bad = make_fixture_interior_();
  sub_bad.K_gI = sub_bad.K_gI + 1i; % complex
  assert_throws_(@() setup_local_schur(sub_bad, opts), 'K_gI complex'); % % ADDED

  fprintf('PASS: test_setup_local_schur\n');
end

% =========================
% Helpers
% =========================

function ensure_project_paths_()
% % ADDED: Make test runnable from any working directory (once tests/ is on path).
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

function sub = make_fixture_interior_()
% % ADDED: Deterministic tiny SPD blocks (nI=2, nG=2), symmetric-consistent coupling.
  K_II = [2, -1; -1, 2];          % SPD
  K_Ig = [1, 0; 2, 1];            % (2x2)
  K_gI = K_Ig.';                  % symmetric-consistent coupling
  K_gg = [3, 1; 1, 2];            % symmetric
  f_I  = [1; -1];
  f_g  = [0.5; 1.5];

  sub = struct('K_II', K_II, 'K_Ig', K_Ig, 'K_gI', K_gI, 'K_gg', K_gg, ...
               'f_I', f_I, 'f_g', f_g);
end

function sub = make_fixture_no_interior_()
% % ADDED: nI=0 case.
  nG   = 2;
  K_II = zeros(0,0);
  K_Ig = zeros(0,nG);
  K_gI = zeros(nG,0);
  K_gg = [4, 1; 1, 2];
  f_I  = zeros(0,1);
  f_g  = [-2; 3];

  sub = struct('K_II', K_II, 'K_Ig', K_Ig, 'K_gI', K_gI, 'K_gg', K_gg, ...
               'f_I', f_I, 'f_g', f_g);
end

function assert_close_(a, b, rtol, atol, msg)
% % ADDED: scale-aware closeness check (inf-norm).
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
% % ADDED: includes a diagnostic print if no error was thrown (Rule 7.3).
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