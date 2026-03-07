% ============================================================
% File: src/common/spectra/plot_spectrum_sorted_overlay.m
% ============================================================
function fig_path = plot_spectrum_sorted_overlay(eig_feti, eig_bddc, case_id, outdir, opts)
%PLOT_SPECTRUM_SORTED_OVERLAY  Plot sorted eigenvalues (FETI-DP vs BDDC) and export to PDF.
%
% Usage:
%   fig_path = plot_spectrum_sorted_overlay(eig_feti, eig_bddc, case_id, outdir);
%   fig_path = plot_spectrum_sorted_overlay(eig_feti, eig_bddc, case_id, outdir, opts);
%
% Inputs:
%   eig_feti : vector of eigenvalues for FETI-DP (will be sorted internally)
%   eig_bddc : vector of eigenvalues for BDDC   (will be sorted internally)
%   case_id  : string used in filename (and as a fallback in title)
%   outdir   : directory where the PDF figure will be saved
%   opts     : optional struct with fields:
%              - yscale      : 'linear' (default) or 'log'
%              - title       : custom title string (highest priority)
%              - mesh_n      : integer mesh size N (for title: N x N)
%              - sub_nx      : integer number of subdomains in x
%              - sub_ny      : integer number of subdomains in y
%              - xlabel      : default 'index i'
%              - ylabel      : default '\lambda_i'
%              - filename    : custom filename (default auto)
%              - fig_visible : 'off' (default) or 'on'
%              - width_px    : default 900
%              - height_px   : default 600
%              - legend_loc  : default 'northeast' (Octave-safe)
%              - paper_fix   : default true (avoid cropping warnings on print)
%
% Output:
%   fig_path : full path to the saved PDF figure.

  % ----------------------------
  % Input validation + defaults
  % ----------------------------
  if nargin < 4 || nargin > 5
    error('plot_spectrum_sorted_overlay:InvalidNargin', ...
          'Expected (eig_feti, eig_bddc, case_id, outdir [, opts]).');
  end

  if nargin < 5 || isempty(opts)
    opts = struct();
  end
  if ~isstruct(opts)
    error('plot_spectrum_sorted_overlay:InvalidOpts', 'opts must be a struct (or []).');
  end

  if ~ischar(case_id) && ~isstring(case_id)
    error('plot_spectrum_sorted_overlay:InvalidCaseId', 'case_id must be a string.');
  end
  case_id = char(case_id);

  if ~ischar(outdir) && ~isstring(outdir)
    error('plot_spectrum_sorted_overlay:InvalidOutdir', 'outdir must be a string.');
  end
  outdir = char(outdir);

  if ~exist(outdir, 'dir')
    mkdir(outdir);
  end

  if ~isfield(opts, 'yscale') || isempty(opts.yscale)
    opts.yscale = 'linear';
  end
  if ~any(strcmp(opts.yscale, {'linear','log'}))
    error('plot_spectrum_sorted_overlay:InvalidYScale', ...
          'opts.yscale must be ''linear'' or ''log''.');
  end

  if ~isfield(opts, 'xlabel') || isempty(opts.xlabel)
    opts.xlabel = 'index i';
  end
  if ~isfield(opts, 'ylabel') || isempty(opts.ylabel)
    opts.ylabel = '\lambda_i';
  end

  if ~isfield(opts, 'fig_visible') || isempty(opts.fig_visible)
    opts.fig_visible = 'off';
  end
  if ~any(strcmp(opts.fig_visible, {'on','off'}))
    error('plot_spectrum_sorted_overlay:InvalidFigVisible', ...
          'opts.fig_visible must be ''on'' or ''off''.');
  end

  if ~isfield(opts, 'width_px') || isempty(opts.width_px)
    opts.width_px = 900;
  end
  if ~isfield(opts, 'height_px') || isempty(opts.height_px)
    opts.height_px = 600;
  end

  if ~isfield(opts, 'legend_loc') || isempty(opts.legend_loc)
    opts.legend_loc = 'northeast';
  end

  if ~isfield(opts, 'paper_fix') || isempty(opts.paper_fix)
    opts.paper_fix = true;
  end
  if ~islogical(opts.paper_fix) || ~isscalar(opts.paper_fix)
    error('plot_spectrum_sorted_overlay:InvalidPaperFix', ...
          'opts.paper_fix must be a logical scalar.');
  end

  if ~isfield(opts, 'filename') || isempty(opts.filename)
    opts.filename = sprintf('fig_spec_sorted_%s.pdf', case_id);
  end

  fig_path = fullfile(outdir, opts.filename);

  % ----------------------------
  % Prepare eigenvalues (sorted)
  % ----------------------------
  eig_feti = eig_feti(:);
  eig_bddc = eig_bddc(:);

  if ~isnumeric(eig_feti) || ~isnumeric(eig_bddc)
    error('plot_spectrum_sorted_overlay:InvalidEigType', ...
          'eig_feti and eig_bddc must be numeric vectors.');
  end

  eig_feti = sort(real(eig_feti), 'ascend');
  eig_bddc = sort(real(eig_bddc), 'ascend');

  n_f = numel(eig_feti);
  n_b = numel(eig_bddc);

  i_f = (1:n_f).';
  i_b = (1:n_b).';

  % ----------------------------
  % Title (priority: opts.title > mesh/subdomain > readable case_id)
  % ----------------------------
  if isfield(opts, 'title') && ~isempty(opts.title)
    title_str = opts.title;
  else
    has_mesh = isfield(opts,'mesh_n') && ~isempty(opts.mesh_n);
    has_sub  = isfield(opts,'sub_nx') && ~isempty(opts.sub_nx) && isfield(opts,'sub_ny') && ~isempty(opts.sub_ny);

    if has_mesh && has_sub
      title_str = sprintf(['Sorted eigenvalues of FETI-DP and BDDC operators ' ...
                           '(Mesh size: %dx%d; Subdomains: %dx%d)'], ...
                           opts.mesh_n, opts.mesh_n, opts.sub_nx, opts.sub_ny);
    else
      readable_id = strrep(case_id, '_', ' ');
      readable_id = strrep(readable_id, 'sub', 'sub ');
      readable_id = strrep(readable_id, 'seed', 'seed ');
      title_str = sprintf('Sorted eigenvalues of FETI-DP and BDDC operators (%s)', readable_id);
    end
  end

  % ----------------------------
  % Create figure
  % ----------------------------
  fig = figure('Visible', opts.fig_visible);
  set(fig, 'Position', [100, 100, opts.width_px, opts.height_px]);

  if opts.paper_fix
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperOrientation', 'landscape');

    pos = get(fig, 'Position');     % [left bottom width_px height_px]
    w_in = pos(3) / 100;
    h_in = pos(4) / 100;

    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 w_in h_in]);
    set(fig, 'PaperSize', [w_in h_in]);
  end

  hold on;
  grid on;

  plot(i_f, eig_feti, '-', 'DisplayName', 'FETI-DP'); % first
  plot(i_b, eig_bddc, '-', 'DisplayName', 'BDDC');    % second

  xlabel(opts.xlabel);
  ylabel(opts.ylabel);
  title(title_str, 'Interpreter', 'none');

  if strcmp(opts.yscale, 'log')
    set(gca, 'YScale', 'log');
  end

  legend('Location', opts.legend_loc);
  xlim([1, max(n_f, n_b)]);

  hold off;

  % ----------------------------
  % Export (robust): PDF first, EPS fallback
  % ----------------------------
  try
    print(fig, fig_path, '-dpdf');
  catch
    [p, name, ~] = fileparts(fig_path);
    eps_path = fullfile(p, [name, '.eps']);

    warning('plot_spectrum_sorted_overlay:PDFExportFailed', ...
            'PDF export failed, falling back to EPS: %s', eps_path);

    print(fig, eps_path, '-depsc2');
  end

  close(fig);
end