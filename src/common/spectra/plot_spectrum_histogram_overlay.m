% ============================================================
% File: src/common/spectra/plot_spectrum_histogram_overlay.m
% ============================================================
function fig_path = plot_spectrum_histogram_overlay(eig_feti, eig_bddc, case_id, outdir, opts)
%PLOT_SPECTRUM_HISTOGRAM_OVERLAY  Compare eigenvalue distributions (FETI-DP vs BDDC) via histogram.
%
% Usage:
%   fig_path = plot_spectrum_histogram_overlay(eig_feti, eig_bddc, case_id, outdir);
%   fig_path = plot_spectrum_histogram_overlay(eig_feti, eig_bddc, case_id, outdir, opts);
%
% Inputs:
%   eig_feti : eigenvalues (vector) for FETI-DP
%   eig_bddc : eigenvalues (vector) for BDDC
%   case_id  : string used in filename (and as a fallback in title)
%   outdir   : directory where the PDF figure will be saved
%   opts     : optional struct with fields:
%              - nbins         : number of bins (default 30)
%              - use_log_bins  : [] (auto) / true / false
%              - log_threshold : ratio threshold for auto log bins (default 1e3)
%              - normalize     : 'count' (default) or 'probability'
%              - title         : custom title string (highest priority)
%              - mesh_n        : for title (Mesh size: NxN)
%              - sub_nx        : for title (Subdomains: sub_nx x sub_ny)
%              - sub_ny        : for title
%              - xlabel        : default '\lambda'
%              - ylabel        : default 'count' or 'probability'
%              - filename      : default 'fig_spec_hist_<case_id>.pdf'
%              - fig_visible   : 'off' (default) or 'on'
%              - width_px      : default 900
%              - height_px     : default 600
%              - legend_loc    : default 'northeast'
%              - paper_fix     : default true (avoid cropping warnings on print)
%
% Output:
%   fig_path : full path to the saved PDF figure.
%
% Notes:
% - We use identical bin edges for both methods (required for fair comparison).
% - Tiny nonpositive eigenvalues (roundoff artifacts) are discarded for the
%   histogram and the discarded counts are reported in the title.

  % ----------------------------
  % Input validation + defaults
  % ----------------------------
  if nargin < 4 || nargin > 5
    error('plot_spectrum_histogram_overlay:InvalidNargin', ...
          'Expected (eig_feti, eig_bddc, case_id, outdir [, opts]).');
  end

  if nargin < 5 || isempty(opts)
    opts = struct();
  end
  if ~isstruct(opts)
    error('plot_spectrum_histogram_overlay:InvalidOpts', 'opts must be a struct (or []).');
  end

  if ~ischar(case_id) && ~isstring(case_id)
    error('plot_spectrum_histogram_overlay:InvalidCaseId', 'case_id must be a string.');
  end
  case_id = char(case_id);

  if ~ischar(outdir) && ~isstring(outdir)
    error('plot_spectrum_histogram_overlay:InvalidOutdir', 'outdir must be a string.');
  end
  outdir = char(outdir);

  if ~exist(outdir, 'dir')
    mkdir(outdir);
  end

  if ~isfield(opts, 'nbins') || isempty(opts.nbins)
    opts.nbins = 30;
  end
  assert(isfinite(opts.nbins) && opts.nbins == floor(opts.nbins) && opts.nbins >= 5, ...
         'opts.nbins must be an integer >= 5.');

  if ~isfield(opts, 'log_threshold') || isempty(opts.log_threshold)
    opts.log_threshold = 1e3;
  end
  assert(isfinite(opts.log_threshold) && opts.log_threshold > 1, ...
         'opts.log_threshold must be > 1.');

  if ~isfield(opts, 'use_log_bins')
    opts.use_log_bins = [];
  end
  if ~(isempty(opts.use_log_bins) || (islogical(opts.use_log_bins) && isscalar(opts.use_log_bins)))
    error('plot_spectrum_histogram_overlay:InvalidUseLogBins', ...
          'opts.use_log_bins must be [] (auto) or a logical scalar.');
  end

  if ~isfield(opts, 'normalize') || isempty(opts.normalize)
    opts.normalize = 'count';
  end
  if ~any(strcmp(opts.normalize, {'count','probability'}))
    error('plot_spectrum_histogram_overlay:InvalidNormalize', ...
          'opts.normalize must be ''count'' or ''probability''.');
  end

  if ~isfield(opts, 'xlabel') || isempty(opts.xlabel)
    opts.xlabel = '\lambda';
  end

  if ~isfield(opts, 'ylabel') || isempty(opts.ylabel)
    if strcmp(opts.normalize, 'probability')
      opts.ylabel = 'probability';
    else
      opts.ylabel = 'count';
    end
  end

  if ~isfield(opts, 'fig_visible') || isempty(opts.fig_visible)
    opts.fig_visible = 'off';
  end
  if ~any(strcmp(opts.fig_visible, {'on','off'}))
    error('plot_spectrum_histogram_overlay:InvalidFigVisible', ...
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
    error('plot_spectrum_histogram_overlay:InvalidPaperFix', ...
          'opts.paper_fix must be a logical scalar.');
  end

  if ~isfield(opts, 'filename') || isempty(opts.filename)
    opts.filename = sprintf('fig_spec_hist_%s.pdf', case_id);
  end

  fig_path = fullfile(outdir, opts.filename);

  % ----------------------------
  % Prepare eigenvalues (filter nonpositive)
  % ----------------------------
  eig_feti = real(eig_feti(:));
  eig_bddc = real(eig_bddc(:));

  if ~isnumeric(eig_feti) || ~isnumeric(eig_bddc)
    error('plot_spectrum_histogram_overlay:InvalidEigType', ...
          'eig_feti and eig_bddc must be numeric vectors.');
  end

  nf_drop = sum(eig_feti <= 0);
  nb_drop = sum(eig_bddc <= 0);

  eig_feti = eig_feti(eig_feti > 0);
  eig_bddc = eig_bddc(eig_bddc > 0);

  if isempty(eig_feti) || isempty(eig_bddc)
    error('plot_spectrum_histogram_overlay:NoPositiveEigenvalues', ...
          'After filtering nonpositive values, at least one eigenvalue set is empty.');
  end

  lam_min = min([eig_feti; eig_bddc]);
  lam_max = max([eig_feti; eig_bddc]);

  if ~(isfinite(lam_min) && isfinite(lam_max) && lam_min > 0 && lam_max > lam_min)
    error('plot_spectrum_histogram_overlay:BadRange', ...
          'Invalid eigenvalue range: min=%.3e, max=%.3e.', lam_min, lam_max);
  end

  ratio = lam_max / lam_min;

  % ----------------------------
  % Decide linear vs log bins
  % ----------------------------
  if isempty(opts.use_log_bins)
    use_log = (ratio >= opts.log_threshold);
  else
    use_log = opts.use_log_bins;
  end

  % ----------------------------
  % Build shared bin edges
  % ----------------------------
  if use_log
    edges = logspace(log10(lam_min), log10(lam_max), opts.nbins + 1);
  else
    edges = linspace(lam_min, lam_max, opts.nbins + 1);
  end

  % ----------------------------
  % Histogram counts (Octave-compatible)
  % ----------------------------
  c_f = histcounts_compat_(eig_feti, edges);
  c_b = histcounts_compat_(eig_bddc, edges);

  centers = 0.5 * (edges(1:end-1) + edges(2:end));

  if strcmp(opts.normalize, 'probability')
    c_f = c_f / max(1, sum(c_f));
    c_b = c_b / max(1, sum(c_b));
  end

  % ----------------------------
  % Title (priority: opts.title > mesh/subdomain > readable case_id)
  % ----------------------------
  title_str = build_title_(opts, case_id);

  if (nf_drop > 0) || (nb_drop > 0)
    title_str = sprintf('%s  [discarded: FETI-DP %d, BDDC %d]', ...
                        title_str, nf_drop, nb_drop);
  end

  % ----------------------------
  % Create figure
  % ----------------------------
  fig = figure('Visible', opts.fig_visible);
  set(fig, 'Position', [100, 100, opts.width_px, opts.height_px]);

  if opts.paper_fix
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperOrientation', 'landscape');

    pos = get(fig, 'Position');   % [left bottom width_px height_px]
    w_in = pos(3) / 100;
    h_in = pos(4) / 100;

    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 w_in h_in]);
    set(fig, 'PaperSize', [w_in h_in]);
  end

  hold on;
  grid on;

  % Grouped bars: each bin center -> [FETI-DP, BDDC]
  Y = [c_f(:), c_b(:)];
  bar(centers, Y, 'grouped');

  xlabel(opts.xlabel);
  ylabel(opts.ylabel);
  title(title_str, 'Interpreter', 'none');

  legend({'FETI-DP', 'BDDC'}, 'Location', opts.legend_loc);

  if use_log
    set(gca, 'XScale', 'log');
  end

  hold off;

  % ----------------------------
  % Export (robust): PDF first, EPS fallback
  % ----------------------------
  try
    print(fig, fig_path, '-dpdf');
  catch
    [p, name, ~] = fileparts(fig_path);
    eps_path = fullfile(p, [name, '.eps']);

    warning('plot_spectrum_histogram_overlay:PDFExportFailed', ...
            'PDF export failed, falling back to EPS: %s', eps_path);

    print(fig, eps_path, '-depsc2');
  end

  close(fig);
end


% ============================================================
% Local helpers
% ============================================================

function counts = histcounts_compat_(x, edges)
%HISTCOUNTS_COMPAT_  Compatibility histogram counts for Octave/MATLAB.
%
% Returns counts for bins [edges(k), edges(k+1)) for k=1..end-1,
% with the last bin including edges(end).

  x = x(:);
  edges = edges(:);

  if exist('histcounts', 'file') == 2
    counts = histcounts(x, edges);
    return;
  end

  % Fallback: histc (older Octave)
  c = histc(x, edges);   %#ok<HISTC>
  % histc returns length(edges) counts; the last count corresponds to x==edges(end)
  % Move that into the last bin and drop the final element.
  if numel(c) ~= numel(edges)
    error('histcounts_compat_:UnexpectedHistcSize', 'histc returned unexpected size.');
  end
  c(end-1) = c(end-1) + c(end);
  counts = c(1:end-1).';
end


function title_str = build_title_(opts, case_id)
%BUILD_TITLE_  Construct a thesis-friendly title with sensible defaults.

  if isfield(opts, 'title') && ~isempty(opts.title)
    title_str = opts.title;
    return;
  end

  has_mesh = isfield(opts,'mesh_n') && ~isempty(opts.mesh_n);
  has_sub  = isfield(opts,'sub_nx') && ~isempty(opts.sub_nx) && ...
             isfield(opts,'sub_ny') && ~isempty(opts.sub_ny);

  if has_mesh && has_sub
    title_str = sprintf('Eigenvalue histogram (Mesh size: %dx%d Subdomains: %dx%d)', ...
                        opts.mesh_n, opts.mesh_n, opts.sub_nx, opts.sub_ny);
    return;
  end

  readable_id = strrep(case_id, '_', ' ');
  readable_id = strrep(readable_id, 'sub', 'sub ');
  readable_id = strrep(readable_id, 'seed', 'seed ');
  title_str = sprintf('Eigenvalue histogram (%s)', readable_id);
end