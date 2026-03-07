function test_bddc_block_integration()
%TEST_BDDC_BLOCK_INTEGRATION  Intermediate/global integration test for the full BDDC block.
%
% Verifies:
%   1) Explicit hat operator A_hat = R^T S R is symmetric and matches applyA_hat().
%   2) Preconditioner applyM_bddc is approximately symmetric (bilinear check) and SPD-ish.
%   3) solve_bddc produces small residual in the explicit A_hat system.
%   4) reconstructed u_free matches monolithic FEM reference Kff\Ff.
%   5) Negative/fault-injection cases throw (overlapping prod_idx; broken SPD in Sdd).
%
% Octave-compatible: plain assertions.

  ensure_project_paths_();

  if exist('rng_deterministic','file') == 2
    rng_deterministic(1);
  end

  f = @(x,y) 2.0 + x - 3*y;  % deterministic RHS (same spirit as DDM integration)

  cases = { ...
    struct('n',  8, 'nSubX', 2, 'nSubY', 2), ...
    struct('n', 12, 'nSubX', 3, 'nSubY', 2) ...
  };

  for ci = 1:numel(cases)
    cfg = cases{ci};
    cfg.f_handle = f;
    cfg.seed = 7 + ci;

    % Keep tolerances tight but realistic for matrix-free vs explicit equality.
    cfg.opts_bddc  = struct('tol_chol', 0);                 % baseline
    cfg.opts_solve = struct('tol', 1e-10, 'maxit', 400);
    cfg.solve = true;

    out = run_bddc_block_case(cfg);
    d = out.diagnostics;

    assert(d.nHat > 0, sprintf('Expected non-empty hat space (case %d).', ci));

    % 1) A_hat symmetry + apply correctness
    assert(d.rel_sym_A < 1e-12, sprintf('A_hat not symmetric enough (case %d): rel=%g', ci, d.rel_sym_A));
    assert(d.rel_applyA < 1e-10, sprintf('applyA_hat mismatch too large (case %d): rel=%g', ci, d.rel_applyA));

    % 2) Preconditioner symmetry/SPD sanity
    assert(d.rel_sym_M < 1e-10, sprintf('applyM_bddc symmetry check failed (case %d): rel=%g', ci, d.rel_sym_M));
    assert(d.M_pos_quad > 0, sprintf('applyM_bddc SPD sanity failed (case %d): x''Mx=%g', ci, d.M_pos_quad));

    % 3) PCG convergence/residual
    assert(out.sol.stats.flag == 0, sprintf('BDDC PCG did not converge (case %d): flag=%d', ci, out.sol.stats.flag));
    assert(d.rel_res_pcg < 1e-8, sprintf('PCG residual too large (case %d): rel=%g', ci, d.rel_res_pcg));

    % (Optional but strong) compare to direct hat solve when small
    if isfinite(d.rel_u_hat_direct)
      assert(d.rel_u_hat_direct < 1e-8, sprintf('u_hat mismatch vs direct solve (case %d): rel=%g', ci, d.rel_u_hat_direct));
    end

    % 4) Reconstruction matches monolithic reference
    assert(d.rel_u_free_ref < 1e-8, sprintf('u_free mismatch vs FEM reference (case %d): rel=%g', ci, d.rel_u_free_ref));

    fprintf('PASS: test_bddc_block_integration case %d (n=%d, %dx%d): iter=%d, rel_u=%.3e\n', ...
            ci, cfg.n, cfg.nSubX, cfg.nSubY, out.sol.stats.iter, d.rel_u_free_ref);
  end

  % ------------------------------------------------------------
  % Negative / fault injection 1: overlapping prod_idx must throw in setup_bddc
  % ------------------------------------------------------------
  cfg = struct('n', 8, 'nSubX', 2, 'nSubY', 2, 'f_handle', f);
  cfg.mutate_data = @mutate_overlap_prod_idx_;
  assert_throws_(@() run_bddc_block_case(cfg), 'overlapping prod_idx should throw');

  % ------------------------------------------------------------
  % Negative / fault injection 2: make Sdd indefinite => chol must fail in setup_bddc
  % ------------------------------------------------------------
  cfg = struct('n', 8, 'nSubX', 2, 'nSubY', 2, 'f_handle', f);
  cfg.mutate_data = @mutate_break_spd_in_Sdd_;
  assert_throws_(@() run_bddc_block_case(cfg), 'non-SPD Sdd should throw');

  fprintf('PASS: test_bddc_block_integration\n');
end

% ========================================================================
% Helpers (same pattern as your DDM integration tests)
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
  if exist('setup_paths','file') == 2
    setup_paths();
  end
  addpath(genpath(fullfile(rootdir, 'tests')));

  did = true;
end

function data = mutate_overlap_prod_idx_(data)
  % Force an overlap between sub(1).prod_idx and sub(2).prod_idx
  if numel(data.sub) < 2
    error('mutate_overlap_prod_idx_: need at least 2 subdomains.');
  end
  a = data.sub(1).prod_idx(:);
  b = data.sub(2).prod_idx(:);
  if isempty(a) || isempty(b)
    error('mutate_overlap_prod_idx_: prod_idx unexpectedly empty.');
  end
  data.sub(2).prod_idx(1) = a(1); % creates overlap
end

function data = mutate_break_spd_in_Sdd_(data)
  % Make Sdd (sub(i).S(idx_d,idx_d)) indefinite by flipping a diagonal entry.
  for i = 1:numel(data.sub)
    if ~isfield(data.sub(i),'S') || isempty(data.sub(i).S), continue; end
    if ~isfield(data.sub(i),'idx_d') || isempty(data.sub(i).idx_d), continue; end
    id = data.sub(i).idx_d(:);
    if isempty(id), continue; end
    j = id(1);
    data.sub(i).S(j,j) = -abs(data.sub(i).S(j,j)) - 1.0;
    return;
  end
  error('mutate_break_spd_in_Sdd_: could not find a subdomain with nonempty idx_d.');
end

function assert_throws_(fh, label)
  threw = false;
  try
    fh();
  catch
    threw = true;
  end
  assert(threw, sprintf('Expected throw: %s', label));
end