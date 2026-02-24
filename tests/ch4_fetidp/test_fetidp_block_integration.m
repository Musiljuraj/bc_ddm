function test_fetidp_block_integration()
%TEST_FETIDP_BLOCK_INTEGRATION  Intermediate/global test for the whole FETI-DP block.
%
% Verifies, for several structured decompositions:
%   1) End-to-end FETI-DP solve + reconstruction matches global reference solve (Kff\Ff).
%   2) Global PDE residual is small (Kff*u_free - Ff).
%   3) Constraint residual on Delta continuity is small (||Bd*wD||).
%   4) Dual equation residual is small (A_lambda*lambda = b).
%   5) A_lambda behaves symmetrically in bilinear form, and energy is positive.
%   6) One negative/fault case throws (invalid partition).
%
% Octave-compatible: plain assertions.

  ensure_project_paths_();

  f = @(x,y) 2.0 + x - 3.0*y;  % deterministic RHS

  cases = { ...
    struct('n', 12, 'nSubX', 3, 'nSubY', 2), ...
    struct('n',  8, 'nSubX', 2, 'nSubY', 2) ...
  };

  for ci = 1:numel(cases)
    cfg = cases{ci};
    cfg.f_handle = f;
    cfg.tol = 1e-10;
    cfg.maxit = 600;
    cfg.setup_opts = struct('tol_rank', 1e-12);
    cfg.seed = 0;

    out = run_fetidp_block_case(cfg);
    d = out.diagnostics;

    % Basic structural sanity
    assert(d.nLambda > 0, sprintf('Expected non-empty lambda space (case %d).', ci));
    assert(d.nDeltaProd >= 0, sprintf('nDeltaProd invalid (case %d).', ci));

    % PCG convergence
    if isstruct(out.stats) && isfield(out.stats,'flag')
      assert(out.stats.flag == 0, sprintf('PCG did not converge (flag=%d) (case %d).', out.stats.flag, ci));
    end

    % Strongest check: compare reconstructed global solution vs direct solve
    assert(d.rel_err_u < 5e-6, sprintf('u_free vs u_ref mismatch too large (case %d): rel=%g', ci, d.rel_err_u));

    % Global residual should be small (tight, but not “direct-solve tight” because PCG tolerance)
    assert(d.rel_res_global < 5e-8, sprintf('Global residual too large (case %d): rel=%g', ci, d.rel_res_global));

    % Constraints on Delta continuity
    assert(d.constraint_norm < 5e-8, sprintf('Constraint residual ||Bd*wD|| too large (case %d): %g', ci, d.constraint_norm));

    % Dual equation residual A_lambda*lambda=b (internal consistency)
    assert(d.rel_res_dual < 5e-7, sprintf('Dual residual too large (case %d): rel=%g', ci, d.rel_res_dual));

    % Symmetry (bilinear) + SPD probe
    assert(d.rel_sym_bilin < 1e-10, sprintf('A_lambda bilinear symmetry probe failed (case %d): rel=%g', ci, d.rel_sym_bilin));
    assert(d.energy_x > 0, sprintf('A_lambda energy probe not positive (case %d): x''Ax=%g', ci, d.energy_x));

    fprintf('PASS: test_fetidp_block_integration case %d (n=%d, %dx%d)\n', ci, cfg.n, cfg.nSubX, cfg.nSubY);
  end

  % ------------------------------------------------------------
  % Negative/fault case: invalid structured partition must throw
  % ------------------------------------------------------------
  cfg_bad = struct('n', 10, 'nSubX', 3, 'nSubY', 2, 'f_handle', f);
  assert_throws_(@() run_fetidp_block_case(cfg_bad), 'invalidPartition should throw');

  fprintf('PASS: test_fetidp_block_integration\n');
end

% ========================================================================
% Helpers (match your established testing style)
% ========================================================================

function ensure_project_paths_()
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

  thisdir  = fileparts(mfilename('fullpath'));
  testsdir = fileparts(thisdir);
  rootdir  = fileparts(testsdir);
  maindir  = fullfile(rootdir, 'main');

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, 'tests')));

  did = true;
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