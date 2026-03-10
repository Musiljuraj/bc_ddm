% ============================================================
% File: main/ch6_run_case.m
% ============================================================
function out = ch6_run_case(cfg)
%CH6_RUN_CASE  Run one Chapter 6 case (FETI-DP vs BDDC) with spectra + PCG stats.
%
% Usage:
%   out = ch6_run_case();   % defaults
%   out = ch6_run_case(struct('n',64,'nSubX',4,'nSubY',4,'seed',1));
%
% Purpose (Chapter 6 workflow):
%   For one small test instance (mesh + subdomain partition), this routine:
%     1) builds the shared Poisson model problem,
%     2) sets up FETI-DP and BDDC on the same instance (FETI-DP first),
%     3) runs PCG to obtain solver diagnostics (iter, relres, resvec),
%     4) computes the full spectra of the preconditioned PCG-seen operators,
%        if the PCG-space dimensions are small enough.
%
% Spectra are computed using:
%   - explicit assembly by basis-vector probing (assemble_from_apply),
%   - Cholesky-based symmetric similarity transform (full_spectrum_precond).
%
% Important control:
%   Full spectra are computed only if n_method <= cfg.nmax, where n_method is:
%     - FETI-DP: nLambda (multiplier space dimension),
%     - BDDC:    nHat = size(R,2) (assembled hat-space dimension).
%
% Output:
%   out : struct with fields:
%         - out.cfg      : used configuration
%         - out.case_id  : deterministic id string for filenames/captions
%         - out.base_meta: basic problem info (no large matrices)
%         - out.fetidp   : { n, stats, spec, spec_skipped, spec_skip_reason }
%         - out.bddc     : { n, stats, spec, spec_skipped, spec_skip_reason }
%
% Notes:
% - By default we do NOT store the full base/data structs to keep .mat small.
% - FETI-DP is always processed first, BDDC second.

  % ----------------------------
  % cfg defaults + validation
  % ----------------------------
  if nargin < 1 || isempty(cfg)
    cfg = struct();
  end
  if ~isstruct(cfg)
    error('ch6_run_case:badCfg', 'cfg must be a struct or empty.');
  end

  % Core case parameters
  if ~isfield(cfg,'n'),      cfg.n = 64;      end
  if ~isfield(cfg,'nSubX'),  cfg.nSubX = 4;   end
  if ~isfield(cfg,'nSubY'),  cfg.nSubY = 4;   end
  if ~isfield(cfg,'seed'),   cfg.seed = 1;    end

  % Solver parameters
  if ~isfield(cfg,'tol'),    cfg.tol = 1e-10; end
  if ~isfield(cfg,'maxit'),  cfg.maxit = 500; end

  % RHS for Poisson (default constant)
  if ~isfield(cfg,'f_handle') || isempty(cfg.f_handle)
    cfg.f_handle = @(x,y) 1.0;
  end

  % Spectra control
  if ~isfield(cfg,'do_spectra'), cfg.do_spectra = true; end
  if ~isfield(cfg,'nmax'),       cfg.nmax = 800;        end

  % Diagnostics / printing
  if ~isfield(cfg,'verbose'), cfg.verbose = false; end

  % Optional: keep base/data in output (default false to keep results small)
  if ~isfield(cfg,'store_base'), cfg.store_base = false; end
  if ~isfield(cfg,'store_data'), cfg.store_data = false; end

  % Optional setup options (kept generic)
  if ~isfield(cfg,'fetidp'), cfg.fetidp = struct(); end
  if ~isfield(cfg.fetidp,'setup_opts'), cfg.fetidp.setup_opts = []; end

  if ~isfield(cfg,'bddc'), cfg.bddc = struct(); end
  if ~isfield(cfg.bddc,'setup_opts'), cfg.bddc.setup_opts = []; end

  % Spectra options passed to full_spectrum_precond
  if ~isfield(cfg,'spectra'), cfg.spectra = struct(); end
  if ~isfield(cfg.spectra,'symmetrize'),     cfg.spectra.symmetrize = true;  end
  if ~isfield(cfg.spectra,'store_matrices'), cfg.spectra.store_matrices = false; end
  if ~isfield(cfg.spectra,'verbose'),        cfg.spectra.verbose = false; end
  if ~isfield(cfg.spectra,'assemble_opts'),  cfg.spectra.assemble_opts = struct(); end

  % Basic numeric validation
  assert(isfinite(cfg.n)      && cfg.n      == floor(cfg.n)      && cfg.n      >= 2, ...
         'cfg.n must be an integer >= 2.');
  assert(isfinite(cfg.nSubX)  && cfg.nSubX  == floor(cfg.nSubX)  && cfg.nSubX  >= 1, ...
         'cfg.nSubX must be an integer >= 1.');
  assert(isfinite(cfg.nSubY)  && cfg.nSubY  == floor(cfg.nSubY)  && cfg.nSubY  >= 1, ...
         'cfg.nSubY must be an integer >= 1.');
  assert(isfinite(cfg.seed)   && cfg.seed   == floor(cfg.seed), ...
         'cfg.seed must be an integer.');
  assert(isfinite(cfg.tol)    && cfg.tol    > 0, ...
         'cfg.tol must be a positive scalar.');
  assert(isfinite(cfg.maxit)  && cfg.maxit  == floor(cfg.maxit) && cfg.maxit >= 0, ...
         'cfg.maxit must be an integer >= 0.');
  assert(isfinite(cfg.nmax)   && cfg.nmax   == floor(cfg.nmax)  && cfg.nmax  >= 1, ...
         'cfg.nmax must be an integer >= 1.');

  if mod(cfg.n, cfg.nSubX) ~= 0 || mod(cfg.n, cfg.nSubY) ~= 0
    error('ch6_run_case:badPartition', ...
          'cfg.n must be divisible by cfg.nSubX and cfg.nSubY (got n=%d, %dx%d).', ...
          cfg.n, cfg.nSubX, cfg.nSubY);
  end

  % ----------------------------
  % Paths + determinism
  % ----------------------------
  if exist('setup_paths','file') ~= 2
    error('ch6_run_case:noSetupPaths', 'setup_paths.m not found on path.');
  end
  setup_paths();

  if exist('rng_deterministic','file') == 2
    rng_deterministic(cfg.seed);
  else
    rng(cfg.seed);
  end

  % ----------------------------
  % Build shared base problem instance
  % ----------------------------
  require_file_('build_problem_data');

  base = build_problem_data(cfg.n, cfg.nSubX, cfg.nSubY, cfg.f_handle);

  base_meta = struct();
  base_meta.n        = cfg.n;
  base_meta.nSubX    = cfg.nSubX;
  base_meta.nSubY    = cfg.nSubY;

  if isfield(base,'p'),    base_meta.nNodes = size(base.p, 1); else, base_meta.nNodes = NaN; end
  if isfield(base,'t'),    base_meta.nElem  = size(base.t, 1); else, base_meta.nElem  = NaN; end
  if isfield(base,'free'), base_meta.nFree  = numel(base.free); else, base_meta.nFree = NaN; end

  % ----------------------------
  % FETI-DP (FIRST)
  % ----------------------------
  require_file_('setup_fetidp');
  require_file_('solve_fetidp');
  require_file_('applyA_lambda');
  require_file_('applyM_lambda');
  require_file_('full_spectrum_precond');

  if cfg.verbose
    fprintf('\n[ch6_run_case] FETI-DP setup (n=%d, %dx%d)\n', cfg.n, cfg.nSubX, cfg.nSubY);
  end

  data_f = setup_call_(@setup_fetidp, base, cfg.fetidp.setup_opts);

  if ~isfield(data_f, 'nLambda')
    error('ch6_run_case:fetidpNoNLambda', 'setup_fetidp output must contain data.nLambda.');
  end
  n_f = data_f.nLambda;

  if cfg.verbose
    fprintf('[ch6_run_case] FETI-DP solve (pcg tol=%.1e, maxit=%d)\n', cfg.tol, cfg.maxit);
  end
  [~, stats_f] = solve_fetidp(data_f, cfg.tol, cfg.maxit);

  spec_f = [];
  spec_f_skipped = false;
  spec_f_skip_reason = '';

  if cfg.do_spectra
    if n_f <= cfg.nmax
      applyA_f    = @(x) applyA_lambda(x, data_f);
      applyMinv_f = @(r) applyM_lambda(r, data_f);

      opts_spec = cfg.spectra;
      spec_f = full_spectrum_precond(applyA_f, applyMinv_f, n_f, opts_spec);
    else
      spec_f_skipped = true;
      spec_f_skip_reason = sprintf('nLambda=%d exceeds nmax=%d (skip full spectrum).', n_f, cfg.nmax);
      spec_f = empty_spec_();
    end
  else
    spec_f_skipped = true;
    spec_f_skip_reason = 'cfg.do_spectra=false (skip full spectrum).';
    spec_f = empty_spec_();
  end

  % ----------------------------
  % BDDC (SECOND)
  % ----------------------------
  require_file_('setup_bddc');
  require_file_('solve_bddc');
  require_file_('applyA_hat');
  require_file_('applyM_bddc');

  if cfg.verbose
    fprintf('\n[ch6_run_case] BDDC setup (n=%d, %dx%d)\n', cfg.n, cfg.nSubX, cfg.nSubY);
  end

  data_b = setup_call_(@setup_bddc, base, cfg.bddc.setup_opts);

  if ~isfield(data_b,'bddc') || ~isfield(data_b.bddc,'R')
    error('ch6_run_case:bddcNoR', 'setup_bddc output must contain data.bddc.R.');
  end
  n_b = size(data_b.bddc.R, 2);

  if cfg.verbose
    fprintf('[ch6_run_case] BDDC solve (pcg tol=%.1e, maxit=%d)\n', cfg.tol, cfg.maxit);
  end
  sol_b = solve_bddc(data_b, struct('tol', cfg.tol, 'maxit', cfg.maxit));

  if ~isfield(sol_b,'stats')
    error('ch6_run_case:bddcNoStats', 'solve_bddc must return a struct with field .stats.');
  end
  stats_b = sol_b.stats;

  spec_b = [];
  spec_b_skipped = false;
  spec_b_skip_reason = '';

  if cfg.do_spectra
    if n_b <= cfg.nmax
      applyA_b    = @(x) applyA_hat(x, data_b);
      applyMinv_b = @(r) applyM_bddc(r, data_b);

      opts_spec = cfg.spectra;
      spec_b = full_spectrum_precond(applyA_b, applyMinv_b, n_b, opts_spec);
    else
      spec_b_skipped = true;
      spec_b_skip_reason = sprintf('nHat=%d exceeds nmax=%d (skip full spectrum).', n_b, cfg.nmax);
      spec_b = empty_spec_();
    end
  else
    spec_b_skipped = true;
    spec_b_skip_reason = 'cfg.do_spectra=false (skip full spectrum).';
    spec_b = empty_spec_();
  end

  % ----------------------------
  % Pack output
  % ----------------------------
  out = struct();
  out.cfg = cfg;
  out.case_id = case_id_(cfg);
  out.base_meta = base_meta;

  if cfg.store_base
    out.base = base;
  end
  if cfg.store_data
    out.data_fetidp = data_f;
    out.data_bddc   = data_b;
  end

  out.fetidp = struct();
  out.fetidp.n               = n_f;
  out.fetidp.stats           = stats_f;
  out.fetidp.spec            = spec_f;
  out.fetidp.spec_skipped    = spec_f_skipped;
  out.fetidp.spec_skip_reason = spec_f_skip_reason;

  out.bddc = struct();
  out.bddc.n                 = n_b;
  out.bddc.stats             = stats_b;
  out.bddc.spec              = spec_b;
  out.bddc.spec_skipped      = spec_b_skipped;
  out.bddc.spec_skip_reason  = spec_b_skip_reason;

  % Optional summary print
  if cfg.verbose
    fprintf('\n=== Chapter 6 case summary: %s ===\n', out.case_id);

    iter_f   = getfield_default_(stats_f, 'iter', NaN);
    relres_f = getfield_default_(stats_f, 'relres', NaN);
    fprintf('FETI-DP: nLambda=%d, iter=%d, relres=%.3e', n_f, iter_f, relres_f);
    if ~out.fetidp.spec_skipped && isfield(spec_f,'kappa')
      fprintf(', kappa=%.3e', spec_f.kappa);
    end
    fprintf('\n');

    iter_b   = getfield_default_(stats_b, 'iter', NaN);
    relres_b = getfield_default_(stats_b, 'relres', NaN);
    fprintf('BDDC   : nHat=%d,   iter=%d, relres=%.3e', n_b, iter_b, relres_b);
    if ~out.bddc.spec_skipped && isfield(spec_b,'kappa')
      fprintf(', kappa=%.3e', spec_b.kappa);
    end
    fprintf('\n\n');
  end
end


% ============================================================
% Local helpers
% ============================================================

function require_file_(fname)
  if exist(fname, 'file') ~= 2
    error('ch6_run_case:missingFile', 'Required function not found on path: %s.m', fname);
  end
end

function data = setup_call_(setup_fun, base, setup_opts)
% Call setup_fun(base [, setup_opts]) while remaining compatible with both
% one-arg and two-arg versions.
  try
    if isempty(setup_opts)
      data = setup_fun(base);
    else
      data = setup_fun(base, setup_opts);
    end
  catch
    % Fallback: some versions may accept only one argument.
    data = setup_fun(base);
  end
end

function sid = case_id_(cfg)
% Deterministic case id for filenames and captions.
  sid = sprintf('n%d_sub%dx%d_seed%d', cfg.n, cfg.nSubX, cfg.nSubY, cfg.seed);
end

function s = empty_spec_()
% Provide a stable empty spec struct when spectra are skipped.
  s = struct();
  s.eigvals   = [];
  s.lmin      = NaN;
  s.lmax      = NaN;
  s.kappa     = NaN;
  s.symmA     = NaN;
  s.symmMinv  = NaN;
  s.symmK     = NaN;
  s.chol_ok   = false;
  s.chol_msg  = '';
end

function v = getfield_default_(s, field, default_val)
% Safe getter for summary printing.
  if isstruct(s) && isfield(s, field)
    v = s.(field);
  else
    v = default_val;
  end
end