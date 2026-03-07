% ============================================================
% File: main/run_ch6_spectral_experiments.m
% ============================================================
function run_ch6_spectral_experiments()
%RUN_CH6_SPECTRAL_EXPERIMENTS  Chapter 6 sequential experiments (Set A + Set B).
%
% This version:
%   - loads the case list from ch6_define_cases.m,
%   - runs all cases,
%   - saves one .mat per case to output/mats/ch6/,
%   - produces (per case):
%       (i)   sorted eigenvalue overlay PDF
%       (ii)  histogram overlay PDF
%       (iii) residual history overlay PDF
%     all saved to output/figures/ch6/,
%   - exports a LaTeX summary table to output/tables/ch6/table_ch6_summary.tex

  % ----------------------------
  % Project root + path setup
  % ----------------------------
  maindir = fileparts(mfilename('fullpath'));
  rootdir = fileparts(maindir);

  addpath(genpath(maindir));
  setup_paths();

  % ----------------------------
  % Output directories
  % ----------------------------
  out_mats   = fullfile(rootdir, 'output', 'mats',    'ch6');
  out_figs   = fullfile(rootdir, 'output', 'figures', 'ch6');
  out_tables = fullfile(rootdir, 'output', 'tables',  'ch6');

  ensure_dir_(out_mats);
  ensure_dir_(out_figs);
  ensure_dir_(out_tables);

  % ----------------------------
  % Global run parameters
  % ----------------------------
  tol        = 1e-10;
  maxit      = 300;

  do_spectra = true;
  nmax       = 1200;      % raise if you want to allow larger full spectra (runtime may grow)

  verbose    = true;

  % ----------------------------
  % Load cases
  % ----------------------------
  cases = ch6_define_cases();
  if isempty(cases)
    error('run_ch6_spectral_experiments:NoCases', 'ch6_define_cases returned no cases.');
  end

  fprintf('\n[run_ch6_spectral_experiments] Running %d case(s)\n\n', numel(cases));

  % Collect outputs for table export
  outs = cell(1, numel(cases));

  % ----------------------------
  % Main loop over cases
  % ----------------------------
  for k = 1:numel(cases)

    cfg = struct();
    cfg.n      = cases(k).n;
    cfg.nSubX  = cases(k).nSubX;
    cfg.nSubY  = cases(k).nSubY;
    cfg.seed   = cases(k).seed;

    cfg.tol    = tol;
    cfg.maxit  = maxit;

    cfg.do_spectra = do_spectra;
    cfg.nmax       = nmax;

    cfg.verbose    = verbose;

    fprintf('============================================================\n');
    fprintf('[run_ch6_spectral_experiments] Case %d/%d: n=%d, sub=%dx%d\n', ...
            k, numel(cases), cfg.n, cfg.nSubX, cfg.nSubY);

    % Run one case
    out = ch6_run_case(cfg);
    outs{k} = out;

    % Save per-case .mat
    mat_name = sprintf('case_%s.mat', out.case_id);
    mat_path = fullfile(out_mats, mat_name);

    save(mat_path, 'out');
    fprintf('[run_ch6_spectral_experiments] saved: %s\n', mat_path);

    % ----------------------------
    % Spectral figures (only if spectra exist for BOTH methods)
    % ----------------------------
    if out.fetidp.spec_skipped || out.bddc.spec_skipped
      fprintf('[run_ch6_spectral_experiments] spectra skipped (no spectral figures).\n');
    else
      % (i) Sorted eigenvalue overlay
      fig_opts = struct();
      fig_opts.yscale      = 'linear';
      fig_opts.fig_visible = 'off';
      fig_opts.mesh_n = cfg.n;
      fig_opts.sub_nx = cfg.nSubX;
      fig_opts.sub_ny = cfg.nSubY;

      fig_path = plot_spectrum_sorted_overlay( ...
        out.fetidp.spec.eigvals, ...
        out.bddc.spec.eigvals, ...
        out.case_id, ...
        out_figs, ...
        fig_opts);

      fprintf('[run_ch6_spectral_experiments] sorted spectrum saved: %s\n', fig_path);

      % (ii) Histogram overlay
      hist_opts = struct();
      hist_opts.fig_visible = 'off';
      hist_opts.mesh_n = cfg.n;
      hist_opts.sub_nx = cfg.nSubX;
      hist_opts.sub_ny = cfg.nSubY;

      hist_path = plot_spectrum_histogram_overlay( ...
        out.fetidp.spec.eigvals, ...
        out.bddc.spec.eigvals, ...
        out.case_id, ...
        out_figs, ...
        hist_opts);

      fprintf('[run_ch6_spectral_experiments] histogram saved: %s\n', hist_path);
    end

    % ----------------------------
    % (iii) Residual history overlay (if resvec exists)
    % ----------------------------
    if has_resvec_(out)
      res_opts = struct();
      res_opts.fig_visible = 'off';
      res_opts.mesh_n = cfg.n;
      res_opts.sub_nx = cfg.nSubX;
      res_opts.sub_ny = cfg.nSubY;

      res_path = plot_residual_history_overlay( ...
        out.fetidp.stats.resvec, ...
        out.bddc.stats.resvec, ...
        out.case_id, ...
        out_figs, ...
        res_opts);

      fprintf('[run_ch6_spectral_experiments] residual history saved: %s\n', res_path);
    else
      fprintf('[run_ch6_spectral_experiments] resvec missing (no residual figure).\n');
    end

    fprintf('============================================================\n\n');
  end

  % ----------------------------
  % Export LaTeX summary table
  % ----------------------------
  table_path = fullfile(out_tables, 'table_ch6_summary.tex');

  table_opts = struct();
  table_opts.include_table_env = true;
  table_opts.booktabs = true;
  table_opts.fontsize_cmd = '\small';

  table_opts.caption = ...
    'Summary of spectral indicators and PCG convergence for FETI-DP and BDDC (Chapter 6).';

  table_opts.label = 'tab:ch6-summary';

  export_ch6_table_tex(outs, table_path, table_opts);

  fprintf('[run_ch6_spectral_experiments] table saved: %s\n', table_path);
  fprintf('[run_ch6_spectral_experiments] DONE\n\n');
end


% ============================================================
% Local helpers
% ============================================================

function ensure_dir_(d)
%ENSURE_DIR_  Create directory if it does not exist.
  if ~exist(d, 'dir')
    mkdir(d);
  end
end

function tf = has_resvec_(out)
%HAS_RESVEC_  True if both methods contain stats.resvec.
  tf = false;

  if ~isstruct(out)
    return;
  end
  if ~isfield(out, 'fetidp') || ~isfield(out, 'bddc')
    return;
  end
  if ~isfield(out.fetidp, 'stats') || ~isfield(out.bddc, 'stats')
    return;
  end
  if ~isstruct(out.fetidp.stats) || ~isstruct(out.bddc.stats)
    return;
  end
  if ~isfield(out.fetidp.stats, 'resvec') || ~isfield(out.bddc.stats, 'resvec')
    return;
  end

  tf = true;
end