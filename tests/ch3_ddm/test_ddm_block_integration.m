function test_ddm_block_integration()
%TEST_DDM_BLOCK_INTEGRATION  Intermediate/global test for the full DDM block (Chapter 3 pipeline).
%
% Verifies, for several structured decompositions:
%   1) Sum of lifted local subdomain matrices reproduces global Dirichlet-reduced system (ddm ordering).
%   2) Assembled-interface Schur complement (from local Schur) matches global Schur computed directly.
%   3) Product interface + primal selection + primal maps are globally consistent.
%   4) A small negative/fault-injection case throws as expected.
%
% Octave-compatible: plain assertions.

  ensure_project_paths_();

  f = @(x,y) 2.0 + x - 3*y;  % deterministic RHS (no manufactured solution needed here)

  cases = { ...
    struct('n', 12, 'nSubX', 3, 'nSubY', 2), ...
    struct('n',  8, 'nSubX', 2, 'nSubY', 2) ...
  };

  for ci = 1:numel(cases)
    cfg = cases{ci};
    cfg.f_handle = f;

    out = run_ddm_block_case(cfg);
    d = out.diagnostics;

    % Must have a nontrivial interface for this DDM block test
    assert(out.iface.nHat > 0, sprintf('Expected non-empty interface for case %d.', ci));

    % ----------------------------
    % 1) Local assembly reproduces global reduced system
    % ----------------------------
    assert(d.rel_Ksum < 1e-12, sprintf('Ksum mismatch too large (case %d): rel=%g', ci, d.rel_Ksum));
    assert(d.rel_Fsum < 1e-12, sprintf('Fsum mismatch too large (case %d): rel=%g', ci, d.rel_Fsum));

    % ----------------------------
    % 2) Schur complement consistency (strongest intermediate check)
    % ----------------------------
    assert(d.rel_sym_S < 1e-12, sprintf('Assembled Schur not symmetric enough (case %d): rel=%g', ci, d.rel_sym_S));

    % Relative tolerances: should be near roundoff; keep conservative but strict.
    assert(d.rel_S < 1e-10, sprintf('Schur mismatch too large (case %d): rel=%g', ci, d.rel_S));
    assert(d.rel_g < 1e-10, sprintf('Reduced RHS mismatch too large (case %d): rel=%g', ci, d.rel_g));

    % ----------------------------
    % 3) Primal selection count (matches your unit-test formula)
    % Dirichlet at x=0 and x=1 => no free nodes there => corners on those lines excluded
    % expected = (nSubX-1)*(nSubY+1)
    % ----------------------------
    expected_nC = (cfg.nSubX - 1) * (cfg.nSubY + 1);
    assert(out.primal.nC == expected_nC, sprintf('Primal count mismatch (case %d): got %d expected %d', ...
           ci, out.primal.nC, expected_nC));

    % ----------------------------
    % 4) Product/primal map global sanity
    % ----------------------------
    assert(out.prod.nProd == sum(arrayfun(@(s) numel(s.gamma_glob), out.sub)), ...
           sprintf('nProd mismatch vs sum |gamma| (case %d).', ci));

    % prod2hat must be consistent with each subdomain's gamma_hat
    for i = 1:out.ddm.Nsub
      pi = out.sub(i).prod_idx(:);
      hi = out.sub(i).gamma_hat(:);
      if isempty(pi), continue; end
      assert(all(out.prod.prod2hat(pi) == hi), sprintf('prod2hat mismatch on sub(%d) (case %d).', i, ci));
    end

    % Each subdomain's idx_c/idx_d should partition its local interface vector
    for i = 1:out.ddm.Nsub
      nG = numel(out.sub(i).gamma_glob);
      ic = out.sub(i).idx_c(:);
      id = out.sub(i).idx_d(:);

      assert(isempty(intersect(ic, id)), sprintf('idx_c and idx_d overlap on sub(%d) (case %d).', i, ci));
      assert(numel(unique([ic; id])) == numel([ic; id]), sprintf('Duplicate indices in idx_c/idx_d on sub(%d) (case %d).', i, ci));
      if nG > 0
        allidx = sort([ic; id]);
        assert(isequal(allidx, (1:nG).'), sprintf('idx_c/idx_d do not cover 1..nG on sub(%d) (case %d).', i, ci));
      else
        assert(isempty(ic) && isempty(id), sprintf('Expected empty idx_c/idx_d for nG=0 on sub(%d) (case %d).', i, ci));
      end
    end

    fprintf('PASS: test_ddm_block_integration case %d (n=%d, %dx%d)\n', ci, cfg.n, cfg.nSubX, cfg.nSubY);
  end

  % ----------------------------
  % Negative / fault-injection: duplicate element index must error in local assembly
  % ----------------------------
  cfg = struct('n', 8, 'nSubX', 2, 'nSubY', 2, 'f_handle', f);
  cfg.mutate_sub = @mutate_duplicate_elems_;
  assert_throws_(@() run_ddm_block_case(cfg), 'duplicate elems should throw');

  fprintf('PASS: test_ddm_block_integration\n');
end

% ========================================================================
% Helpers
% ========================================================================

function ensure_project_paths_()
  % Same “run anywhere” style used in your unit tests.
  persistent did;
  if ~isempty(did) && did
    return;
  end

  if exist('setup_paths','file') == 2
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

  thisdir  = fileparts(mfilename('fullpath'));   % .../tests/ch3_ddm
  testsdir = fileparts(thisdir);                 % .../tests
  rootdir  = fileparts(testsdir);                % project root
  maindir  = fullfile(rootdir, 'main');

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, 'tests')));

  did = true;
end

function [sub, ddm] = mutate_duplicate_elems_(sub, ddm)
  %#ok<INUSD>
  % Make a guaranteed duplicate in sub(1).elems (if possible)
  if ~isempty(sub) && isfield(sub(1),'elems') && numel(sub(1).elems) >= 1
    sub(1).elems(end+1,1) = sub(1).elems(1);
  end
end

function assert_throws_(fh, label)
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  assert(threw, sprintf('Expected an error, but none was thrown (%s).', label));
end