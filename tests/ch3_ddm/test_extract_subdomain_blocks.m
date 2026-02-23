function test_extract_subdomain_blocks()
%TEST_EXTRACT_SUBDOMAIN_BLOCKS  Deterministic unit tests for extract_subdomain_blocks.
%
% Run (from repo root):
%   addpath(genpath('tests')); clear functions; rehash; test_extract_subdomain_blocks
%
% % ADDED: path-robust setup so the test can run from any working directory
  ensure_project_paths_();

  % -----------------------------
  % % CHANGED: tiny deterministic fixtures (no FEM/DDM pipeline dependencies)
  % -----------------------------

  % Case 1: single subdomain, dense K, nontrivial I/G split
  sub1 = struct();
  sub1.K = [ 4,  1, -2;
             1,  3,  0;
            -2,  0,  5 ];
  sub1.f = [10; 20; 30];
  sub1.dofs_I = [1; 3];
  sub1.dofs_G = 2;

  out = extract_subdomain_blocks(sub1);

  % % ADDED: field presence + size/type invariants
  assert(isstruct(out), 'Expected struct output.');
  assert(isfield(out, 'K_II') && isfield(out, 'K_Ig') && isfield(out, 'K_gI') && isfield(out, 'K_gg'), ...
         'Missing expected K_* fields.');
  assert(isfield(out, 'f_I')  && isfield(out, 'f_g'), 'Missing expected f_* fields.');
  assert(isfield(out, 'nI')   && isfield(out, 'nG'),  'Missing expected nI/nG fields.');

  I = sub1.dofs_I(:); G = sub1.dofs_G(:);
  assert(out.nI == numel(I), 'nI does not match numel(dofs_I).');
  assert(out.nG == numel(G), 'nG does not match numel(dofs_G).');

  assert(isequal(size(out.K_II), [numel(I), numel(I)]), 'K_II size mismatch.');
  assert(isequal(size(out.K_Ig), [numel(I), numel(G)]), 'K_Ig size mismatch.');
  assert(isequal(size(out.K_gI), [numel(G), numel(I)]), 'K_gI size mismatch.');
  assert(isequal(size(out.K_gg), [numel(G), numel(G)]), 'K_gg size mismatch.');

  assert(isequal(size(out.f_I), [numel(I), 1]), 'f_I size mismatch.');
  assert(isequal(size(out.f_g), [numel(G), 1]), 'f_g size mismatch.');

  % % ADDED: exact reconstruction check (slicing should match exactly)
  perm = [I; G];
  Kp   = sub1.K(perm, perm);
  Krec = [out.K_II, out.K_Ig;
          out.K_gI, out.K_gg];
  assert(isequal(Kp, Krec), 'Block reconstruction mismatch (Case 1).');

  fp   = sub1.f(perm);
  frec = [out.f_I; out.f_g];
  assert(isequal(fp, frec), 'RHS reconstruction mismatch (Case 1).');

  % -----------------------------
  % % FIXED: Octave struct-array construction (avoid "incompatible fields" error)
  % -----------------------------
  tmpl = struct('K', [], 'f', [], 'dofs_I', [], 'dofs_G', []);
  sub  = repmat(tmpl, 1, 2);  % % FIXED
  sub(1) = sub1;

  % Case 2: sparse K, different split, multiple subdomains
  sub(2) = struct( ...
    'K', sparse([ 2, -1,  0,  0;
                 -1,  2, -1,  0;
                  0, -1,  2, -1;
                  0,  0, -1,  2 ]), ...
    'f', [1; 2; 3; 4], ...
    'dofs_I', [2; 3], ...
    'dofs_G', [1; 4] );

  out2 = extract_subdomain_blocks(sub);

  % % ADDED: per-subdomain reconstruction + RHS checks
  for i = 1:numel(sub)
    I = sub(i).dofs_I(:); G = sub(i).dofs_G(:);
    perm = [I; G];

    Kp   = sub(i).K(perm, perm);
    Krec = [out2(i).K_II, out2(i).K_Ig;
            out2(i).K_gI, out2(i).K_gg];
    assert(isequal(Kp, Krec), sprintf('Block reconstruction mismatch (Case 2 sub(%d)).', i));

    fp   = sub(i).f(perm);
    frec = [out2(i).f_I; out2(i).f_g];
    assert(isequal(fp, frec), sprintf('RHS reconstruction mismatch (Case 2 sub(%d)).', i));
  end

  % % ADDED: for sparse input, blocks should remain sparse (Octave indexing preserves sparsity)
  assert(issparse(out2(2).K_II) && issparse(out2(2).K_Ig) && issparse(out2(2).K_gI) && issparse(out2(2).K_gg), ...
         'Expected sparse block outputs for sparse input K.');

  % -----------------------------
  % % ADDED: empty interior / empty interface edge cases
  % -----------------------------
  subE = struct();
  subE.K = [1, 2;
            3, 4];
  subE.f = [5; 6];

  % Empty interior, all interface
  subE.dofs_I = [];
  subE.dofs_G = [1; 2];
  outE = extract_subdomain_blocks(subE);
  assert(outE.nI == 0 && outE.nG == 2, 'Empty-interior nI/nG mismatch.');
  assert(isequal(size(outE.K_II), [0, 0]), 'Empty-interior K_II size mismatch.');
  assert(isequal(size(outE.K_Ig), [0, 2]), 'Empty-interior K_Ig size mismatch.');
  assert(isequal(size(outE.K_gI), [2, 0]), 'Empty-interior K_gI size mismatch.');
  assert(isequal(size(outE.K_gg), [2, 2]), 'Empty-interior K_gg size mismatch.');
  assert(isequal(outE.f_I, zeros(0,1)), 'Empty-interior f_I mismatch.');
  assert(isequal(outE.f_g, subE.f), 'Empty-interior f_g mismatch.');

  % Empty interface, all interior
  subE.dofs_I = [1; 2];
  subE.dofs_G = [];
  outE2 = extract_subdomain_blocks(subE);
  assert(outE2.nI == 2 && outE2.nG == 0, 'Empty-interface nI/nG mismatch.');
  assert(isequal(size(outE2.K_II), [2, 2]), 'Empty-interface K_II size mismatch.');
  assert(isequal(size(outE2.K_Ig), [2, 0]), 'Empty-interface K_Ig size mismatch.');
  assert(isequal(size(outE2.K_gI), [0, 2]), 'Empty-interface K_gI size mismatch.');
  assert(isequal(size(outE2.K_gg), [0, 0]), 'Empty-interface K_gg size mismatch.');
  assert(isequal(outE2.f_I, subE.f), 'Empty-interface f_I mismatch.');
  assert(isequal(outE2.f_g, zeros(0,1)), 'Empty-interface f_g mismatch.');

  % -----------------------------
  % % ADDED: signature tests (wrong arg count / too many outputs must error)
  % -----------------------------
  assert_throws_(@() extract_subdomain_blocks(), ...
                'Expected an error for missing input argument.');
  assert_throws_(@() extract_subdomain_blocks(sub1, 123), ...
                'Expected an error for too many input arguments.');
  assert_throws_(@() call_two_outputs_(sub1), ...
                'Expected an error for too many output arguments.');

  % -----------------------------
  % % ADDED: missing required fields must error
  % -----------------------------
  base = struct('K', eye(2), 'f', [1; 2], 'dofs_I', 1, 'dofs_G', 2);

  tmp = base; tmp = rmfield(tmp, 'K');
  assert_throws_(@() extract_subdomain_blocks(tmp), 'Expected an error when field K is missing.');

  tmp = base; tmp = rmfield(tmp, 'f');
  assert_throws_(@() extract_subdomain_blocks(tmp), 'Expected an error when field f is missing.');

  tmp = base; tmp = rmfield(tmp, 'dofs_I');
  assert_throws_(@() extract_subdomain_blocks(tmp), 'Expected an error when field dofs_I is missing.');

  tmp = base; tmp = rmfield(tmp, 'dofs_G');
  assert_throws_(@() extract_subdomain_blocks(tmp), 'Expected an error when field dofs_G is missing.');

  % -----------------------------
  % % ADDED: invalid/ambiguous index rejection tests (drives function hardening)
  % Rule 7.3 diagnostic prints which bad-case did not throw.
  % -----------------------------
  Kbig = speye(100);      % large enough that ASCII indices like 'A' (65) are in-range
  fbig = zeros(100, 1);

  bad = {};

  % char indices (ASCII indexing)
  s = struct('K', Kbig, 'f', fbig, 'dofs_I', 'A', 'dofs_G', 1);
  bad{end+1} = struct('name', 'dofs_I is char (ASCII indexing)', 'sub', s);

  % logical indices (mask indexing)
  s = struct('K', Kbig, 'f', fbig, 'dofs_I', true(100,1), 'dofs_G', [1; 2]);
  bad{end+1} = struct('name', 'dofs_I is logical (mask indexing)', 'sub', s);

  % duplicate indices
  s = struct('K', Kbig, 'f', fbig, 'dofs_I', [1; 1; 2], 'dofs_G', 3);
  bad{end+1} = struct('name', 'dofs_I contains duplicates', 'sub', s);

  % overlap between interior and interface
  s = struct('K', Kbig, 'f', fbig, 'dofs_I', [1; 2; 3], 'dofs_G', [3; 4]);
  bad{end+1} = struct('name', 'dofs_I and dofs_G overlap', 'sub', s);

  for k = 1:numel(bad)
    threw = false;
    try
      extract_subdomain_blocks(bad{k}.sub);
    catch
      threw = true;
    end
    if ~threw
      fprintf('DIAGNOSTIC: bad-case %d did not throw: %s\n', k, bad{k}.name);
    end
    assert(threw, sprintf('Expected an error but none was thrown for bad-case %d: %s', k, bad{k}.name));
  end

  fprintf('PASS: test_extract_subdomain_blocks\n');
end

% -----------------------------
% % ADDED: helpers
% -----------------------------
function assert_throws_(fh, msg)
  if nargin < 2
    msg = 'Expected an error but none was thrown.';
  end
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  assert(threw, msg);
end

function call_two_outputs_(s)
  % % ADDED: helper to force "too many outputs" path
  [a, b] = extract_subdomain_blocks(s); %#ok<NASGU>
end

function ensure_project_paths_()
%ENSURE_PROJECT_PATHS_  Make the test runnable from any working directory.
%
% Important:
% - Octave must still be able to find this test function initially (tests/ on path once per session).
% - After that, this helper makes the rest deterministic.

  % % ADDED: if setup_paths is already on path, use it and infer root from its location
  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    maindir = fileparts(sp);
    rootdir = fileparts(maindir);
    testsdir = fullfile(rootdir, 'tests');
    if exist(testsdir, 'dir') == 7
      addpath(genpath(testsdir));
    end
    return;
  end

  % % ADDED: locate repo root relative to this test file
  thisdir  = fileparts(mfilename('fullpath'));
  testsdir = fileparts(thisdir);
  rootdir  = fileparts(testsdir);
  maindir  = fullfile(rootdir, 'main');

  if exist(maindir, 'dir') ~= 7
    error('ensure_project_paths_: cannot locate bc_ddm root (expected main/ at %s).', maindir);
  end

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, 'tests')));
end