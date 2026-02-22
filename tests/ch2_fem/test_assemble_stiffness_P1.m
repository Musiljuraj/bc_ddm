function test_assemble_stiffness_P1()
%TEST_ASSEMBLE_STIFFNESS_P1 Strong regression tests for global stiffness assembly.
%
% FIXED: tests are now fully non-interactive (disable recursive rmdir prompts).
% CHANGED/ADDED: broadened coverage (arg-count, invalid inputs, multiple n, invariances).
% ADDED: mapping test that catches Ke vs Ke' assembly using a non-symmetric stub.

  % ADDED: ensure rmdir(...) never prompts "yes or no" during tests
  old_confirm = maybe_disable_recursive_rmdir_confirm();  % ADDED
  cleanup_obj = onCleanup(@() maybe_restore_recursive_rmdir_confirm(old_confirm)); % ADDED

  % Keep tests runnable standalone (harmless if caller already set paths).
  if exist('setup_paths', 'file') == 2
    setup_paths();
  end

  % -----------------------------
  % 0) Signature / arg-count
  % -----------------------------
  % ADDED
  assert_throws(@() assemble_stiffness_P1(), 'nargin=0');
  assert_throws(@() assemble_stiffness_P1([0 0]), 'nargin=1');
  assert_throws(@() assemble_stiffness_P1([0 0], [1 1 1], 123), 'nargin=3');

  % -----------------------------
  % 1) Basic input validation
  % -----------------------------
  % ADDED
  p3 = [0 0; 1 0; 0 1];
  t1 = [1 2 3];

  assert_throws(@() assemble_stiffness_P1(zeros(3,3), t1), 'p not Nx2');
  assert_throws(@() assemble_stiffness_P1(p3, ones(1,4)), 't not Nex3');

  assert_throws(@() assemble_stiffness_P1(p3, [0 2 3]), 't has 0 index');
  assert_throws(@() assemble_stiffness_P1(p3, [1 2 4]), 't has out-of-range index');
  assert_throws(@() assemble_stiffness_P1(p3, [1.1 2 3]), 't has non-integer index');

  % -----------------------------
  % 2) Empty triangulation behavior
  % -----------------------------
  % ADDED
  K0 = assemble_stiffness_P1(zeros(0,2), zeros(0,3));
  assert(issparse(K0), 'Empty case: K should be sparse.');
  assert(all(size(K0) == [0 0]), 'Empty case: K should be 0x0.');

  % -----------------------------
  % 3) Local-to-global mapping correctness (Ke vs Ke' bug catcher)
  % -----------------------------
  % ADDED: This test intentionally uses a non-symmetric Ke stub.
  check_mapping_with_nonsymmetric_stub();  % ADDED

  % -----------------------------
  % 4) Invariants on small structured meshes
  % -----------------------------
  n_list = [1 2 4]; % CHANGED: test multiple sizes

  for ni = 1:numel(n_list)
    n = n_list(ni);

    [p, t, bnd] = mesh_unit_square_P1(n); %#ok<NASGU>
    K = assemble_stiffness_P1(p, t);

    N = size(p,1);

    % size + sparse type
    assert(all(size(K) == [N, N]), sprintf('K has wrong size for n=%d.', n));
    assert(issparse(K), sprintf('K should be sparse for n=%d.', n));

    % ADDED: finite entries
    v = nonzeros(K);
    assert(all(isfinite(v)), sprintf('K has NaN/Inf entries for n=%d.', n));

    % symmetry (scale-aware)
    rel_sym = norm(K - K.', 'fro') / max(1, norm(K, 'fro'));
    assert(rel_sym < 1e-12, sprintf('K not symmetric for n=%d (rel=%g).', n, rel_sym));

    % constant nullspace / row-sum (scale-aware)
    one = ones(N,1);
    r1  = K * one;
    rs  = full(sum(K,2));

    tol_null = scaled_tol_vec(K, N);
    assert(norm(r1, 2) <= tol_null, ...
           sprintf('K*ones not ~0 for n=%d (||r||=%g, tol=%g).', n, norm(r1,2), tol_null));
    assert(norm(rs, 2) <= tol_null, ...
           sprintf('Row sums not ~0 for n=%d (||s||=%g, tol=%g).', n, norm(rs,2), tol_null));

    % diagonal should be nonnegative (allow tiny roundoff)
    tol_d = 10 * eps * max(1, norm(K,1));
    assert(all(diag(K) >= -tol_d), sprintf('Negative diagonal entries for n=%d.', n));

    % ADDED: PSD checks on mean-zero vectors (avoid nullspace)
    rand('state', 1000 + n);  %#ok<RAND>
    randn('state', 2000 + n); %#ok<RANDN>

    for j = 1:3
      x = randn(N,1);
      x = x - mean(x);

      q = x' * (K * x);
      tol_psd = 1e2 * eps * max(1, norm(K,1)) * (norm(x,2)^2);

      assert(q >= -tol_psd, ...
             sprintf('PSD violated (tiny negative energy) for n=%d: q=%g, tol=%g.', n, q, tol_psd));

      % energy invariant under adding constant vector
      c  = 1.2345;
      q2 = (x + c*one)' * (K * (x + c*one));
      tol_q = 1e2 * eps * max(1, abs(q)) + tol_psd;
      assert(abs(q2 - q) <= tol_q, ...
             sprintf('Energy not invariant under constant shift for n=%d.', n));
    end

    % invariance under local vertex permutation
    t_perm = t(:, [2 3 1]);
    Kp = assemble_stiffness_P1(p, t_perm);

    rel_perm = norm(K - Kp, 'fro') / max(1, norm(K, 'fro'));
    assert(rel_perm < 1e-12, sprintf('Not invariant to local vertex permutation for n=%d.', n));

    % invariance under global node renumbering: K_new == P' * K * P
    perm = randperm(N);
    invperm = zeros(N,1); invperm(perm) = 1:N;

    p2 = p(perm, :);
    t2 = invperm(t);

    K2 = assemble_stiffness_P1(p2, t2);

    P = sparse(perm, 1:N, 1, N, N);
    K2_expected = P' * K * P;

    rel_glob = norm(K2 - K2_expected, 'fro') / max(1, norm(K, 'fro'));
    assert(rel_glob < 1e-12, sprintf('Not consistent under global renumbering for n=%d.', n));
  end

  fprintf('PASS: test_assemble_stiffness_P1\n');
end

% -------------------------------------------------------------------------
% Helper: disable/restore recursive rmdir confirmation (non-interactive tests)
% -------------------------------------------------------------------------
function old = maybe_disable_recursive_rmdir_confirm()
  % ADDED
  old = [];
  if exist('confirm_recursive_rmdir', 'file') == 2 || exist('confirm_recursive_rmdir', 'builtin') == 5
    try
      old = confirm_recursive_rmdir();
      confirm_recursive_rmdir(false);
    catch
      old = [];
    end
  end
end

function maybe_restore_recursive_rmdir_confirm(old)
  % ADDED
  if ~isempty(old)
    if exist('confirm_recursive_rmdir', 'file') == 2 || exist('confirm_recursive_rmdir', 'builtin') == 5
      try
        confirm_recursive_rmdir(old);
      catch
        % ignore
      end
    end
  end
end

% -------------------------------------------------------------------------
% Helper: assert a call throws (Octave-portable; no message substring checks)
% -------------------------------------------------------------------------
function assert_throws(fh, label)
  % ADDED
  if nargin < 2
    label = '(no label)';
  end
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  assert(threw, sprintf('Expected an error but none was thrown: %s', label));
end

% -------------------------------------------------------------------------
% Helper: scale-aware tolerance for "near-zero" vector invariants
% -------------------------------------------------------------------------
function tol = scaled_tol_vec(K, N)
  % ADDED
  s = max(1, norm(K, 1));
  tol = max(1e-12 * s * sqrt(max(1,N)), 50 * eps * s * sqrt(max(1,N)));
end

% -------------------------------------------------------------------------
% Helper: mapping test that detects assembling Ke' instead of Ke
% -------------------------------------------------------------------------
function check_mapping_with_nonsymmetric_stub()
  % ADDED
  tmpdir = tempname();
  mkdir(tmpdir);

  stubfile = fullfile(tmpdir, 'triP1_stiffness.m');
  fid = fopen(stubfile, 'w');
  if fid < 0
    error('Failed to create temporary triP1_stiffness stub.');
  end

  % Non-symmetric on purpose (so Ke vs Ke' differs).
  fprintf(fid, 'function Ke = triP1_stiffness(xy)\n');
  fprintf(fid, '%% Temporary stub for assembly mapping test.\n');
  fprintf(fid, '%% xy is unused except for signature compatibility.\n');
  fprintf(fid, 'Ke = [1 2 3; 4 5 6; 7 8 9];\n');
  fprintf(fid, 'end\n');
  fclose(fid);

  oldpath = path();
  try
    addpath(tmpdir, '-begin');

    clear triP1_stiffness;
    clear assemble_stiffness_P1;

    p = [0 0; 1 0; 0 1];
    t = [1 2 3];

    Ke = triP1_stiffness(p(t,:));
    K  = assemble_stiffness_P1(p, t);

    A = full(K);
    err = norm(A - Ke, 'fro');
    assert(err == 0, sprintf(['Assembly mapping error: expected K == Ke for single element.\n' ...
                              'This typically indicates swapped (row,col) triplets (Ke'' assembled).\n' ...
                              '||K-Ke||_F = %g'], err));
  catch err
    % FIXED: explicit cleanup without prompting
    path(oldpath);
    if exist(stubfile, 'file') == 2
      delete(stubfile);
    end
    if exist(tmpdir, 'dir') == 7
      rmdir(tmpdir, 's');
    end
    rethrow(err);
  end

  % FIXED: explicit cleanup without prompting
  path(oldpath);
  if exist(stubfile, 'file') == 2
    delete(stubfile);
  end
  if exist(tmpdir, 'dir') == 7
    rmdir(tmpdir, 's');
  end
  clear triP1_stiffness;
  clear assemble_stiffness_P1;
end