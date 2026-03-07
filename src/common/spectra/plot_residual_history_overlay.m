% ============================================================
% File: src/common/spectra/plot_residual_history_overlay.m
% ============================================================
function fig_path = plot_residual_history_overlay( ...
    resvec_feti, resvec_bddc, case_id, outdir, opts)
%PLOT_RESIDUAL_HISTORY_OVERLAY  Plot PCG residual histories and export to PDF.
%
% Usage:
%   fig_path = plot_residual_history_overlay(res_f, res_b, case_id, outdir);
%   fig_path = plot_residual_history_overlay(res_f, res_b, case_id, outdir, opts);
%
% Inputs:
%   resvec_feti : residual history for FETI-DP (vector)
%   resvec_bddc : residual history for BDDC (vector)
%   case_id     : string used in filename (and fallback title)
%   outdir      : directory for output figure
%   opts        : optional struct:
%                 - title       : custom title (highest priority)
%                 - mesh_n      : for title (Mesh size: NxN)
%                 - sub_nx      : for title (Subdomains: sub_nx x sub_ny)
%                 - sub_ny      : for title
%                 - xlabel      : default 'iteration k'
%                 - ylabel      : default '||r_k|| / ||r_0||'
%                 - filename    : default 'fig_residual_<case_id>.pdf'
%                 - fig_visible : 'off' (default) or 'on'
%                 - width_px    : default 900
%                 - height_px   : default 600
%                 - legend_loc  : default 'northeast'
%                 - paper_fix   : default true
%                 - normalize   : true (default), plot ||r_k||/||r_0||
%
% Output:
%   fig_path : full path to saved PDF figure.

  % ----------------------------
  % Input validation + defaults
  % ----------------------------
  if nargin < 4 || nargin > 5
    error('plot_residual_history_overlay:InvalidNargin', ...
          'Expected (resvec_feti, resvec_bddc, case_id, outdir [, opts]).');
  end

  if nargin < 5 || isempty(opts)
    opts = struct();
  end
  if ~isstruct(opts)
    error('plot_residual_history_overlay:InvalidOpts', 'opts must be a struct (or []).');
  end

  if ~ischar(case_id) && ~isstring(case_id)
    error('plot_residual_history_overlay:InvalidCaseId', 'case_id must be a string.');
  end
  case_id = char(case_id);

  if ~ischar(outdir) && ~isstring(outdir)
    error('plot_residual_history_overlay:InvalidOutdir', 'outdir must be a string.');
  end
  outdir = char(outdir);

  if ~exist(outdir, 'dir')
    mkdir(outdir);
  end

  if ~isfield(opts, 'xlabel') || isempty(opts.xlabel)
    opts.xlabel = 'iteration k';
  end

  if ~isfield(opts, 'ylabel') || isempty(opts.ylabel)
    opts.ylabel = '||r_k|| / ||r_0||';
  end

  if ~isfield(opts, 'fig_visible') || isempty(opts.fig_visible)
    opts.fig_visible = 'off';
  end
  if ~any(strcmp(opts.fig_visible, {'on','off'}))
    error('plot_residual_history_overlay:InvalidFigVisible', ...
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
    error('plot_residual_history_overlay:InvalidPaperFix', ...
          'opts.paper_fix must be a logical scalar.');
  end

  if ~isfield(opts, 'normalize') || isempty(opts.normalize)
    opts.normalize = true;
  end
  if ~islogical(opts.normalize) || ~isscalar(opts.normalize)
    error('plot_residual_history_overlay:InvalidNormalize', ...
          'opts.normalize must be a logical scalar.');
  end

  if ~isfield(opts, 'filename') || isempty(opts.filename)
    opts.filename = sprintf('fig_residual_%s.pdf', case_id);
  end

  fig_path = fullfile(outdir, opts.filename);

  % ----------------------------
  % Prepare residual histories
  % ----------------------------
  resvec_feti = real(resvec_feti(:));
  resvec_bddc = real(resvec_bddc(:));

  if ~isnumeric(resvec_feti) || ~isnumeric(resvec_bddc)
    error('plot_residual_history_overlay:InvalidResvecType', ...
          'resvec_feti and resvec_bddc must be numeric vectors.');
  end

  if isempty(resvec_feti) || isempty(resvec_bddc)
    error('plot_residual_history_overlay:EmptyResvec', ...
          'Residual history vectors must be non-empty.');
  end

  if opts.normalize
    if resvec_feti(1) ~= 0
      resvec_feti = resvec_feti / resvec_feti(1);
    end
    if resvec_bddc(1) ~= 0
      resvec_bddc = resvec_bddc / resvec_bddc(1);
    end
  end

  k_f = (0:(numel(resvec_feti)-1)).';
  k_b = (0:(numel(resvec_bddc)-1)).';

  % ----------------------------
  % Title construction
  % ----------------------------
  if isfield(opts, 'title') && ~isempty(opts.title)
    title_str = opts.title;
  else
    has_mesh = isfield(opts,'mesh_n') && ~isempty(opts.mesh_n);
    has_sub  = isfield(opts,'sub_nx') && ~isempty(opts.sub_nx) && ...
               isfield(opts,'sub_ny') && ~isempty(opts.sub_ny);

    if has_mesh && has_sub
      title_str = sprintf( ...
        'PCG residual histories (Mesh size: %dx%d Subdomains: %dx%d)', ...
        opts.mesh_n, opts.mesh_n, opts.sub_nx, opts.sub_ny);
    else
      readable_id = strrep(case_id, '_', ' ');
      readable_id = strrep(readable_id, 'sub', 'sub ');
      readable_id = strrep(readable_id, 'seed', 'seed ');
      title_str = sprintf('PCG residual histories (%s)', readable_id);
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

  % semilogy is standard for residual histories
  semilogy(k_f, resvec_feti, '-', 'DisplayName', 'FETI-DP');
  semilogy(k_b, resvec_bddc, '-', 'DisplayName', 'BDDC');

  xlabel(opts.xlabel);
  ylabel(opts.ylabel);
  title(title_str, 'Interpreter', 'none');

  legend('Location', opts.legend_loc);
  xlim([0, max(k_f(end), k_b(end))]);

  hold off;

  % ----------------------------
  % Export (robust): PDF first, EPS fallback
  % ----------------------------
  try
    print(fig, fig_path, '-dpdf');
  catch
    [p, name, ~] = fileparts(fig_path);
    eps_path = fullfile(p, [name, '.eps']);

    warning('plot_residual_history_overlay:PDFExportFailed', ...
            'PDF export failed, falling back to EPS: %s', eps_path);

    print(fig, eps_path, '-depsc2');
  end

  close(fig);
end