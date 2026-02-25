% ============================================================
% File: src/bddc/solve/solve_bddc.m
% ============================================================
function out = solve_bddc(data, opts)
%SOLVE_BDDC  Solve the BDDC hat-space interface problem by PCG.
%
%   out = solve_bddc(data, opts)
%
% Solves:
%   A_hat u_hat = f_hat,
%   A_hat = R^T S R,  f_hat = R^T g,
% with baseline BDDC preconditioner.
%
% Inputs:
%   data : struct from build_problem_data(...) extended by setup_bddc(...)
%   opts : optional struct
%          opts.tol   (default 1e-8)
%          opts.maxit (default 500)
%
% Output struct out:
%   out.u_hat : hat-space solution
%   out.stats : pcg_wrap stats
%   out.w_prod: product-space interface trace R*u_hat

  %% HARDENING CHANGE: basic arg-count validation
  if nargin < 1 || nargin > 2
    error('solve_bddc:bad_nargin', 'solve_bddc expects 1 or 2 input arguments.');
  end

  %% HARDENING CHANGE: validate data/opts types and defaults
  if ~isstruct(data)
    error('solve_bddc:bad_data', 'data must be a struct (from setup_bddc).');
  end

  if nargin < 2 || isempty(opts)
    opts = struct();
  end
  if ~isstruct(opts)
    error('solve_bddc:bad_opts', 'opts must be a struct or empty.');
  end
  if ~isfield(opts,'tol');   opts.tol = 1e-8; end
  if ~isfield(opts,'maxit'); opts.maxit = 500; end

  %% HARDENING CHANGE: validate opts fields
  if ~isscalar(opts.tol) || ~isnumeric(opts.tol) || ~(opts.tol > 0) || ~isfinite(opts.tol)
    error('solve_bddc:bad_tol', 'opts.tol must be a positive finite scalar.');
  end
  if ~isscalar(opts.maxit) || ~isnumeric(opts.maxit) || ~(opts.maxit >= 0) || ~isfinite(opts.maxit)
    error('solve_bddc:bad_maxit', 'opts.maxit must be a finite scalar >= 0.');
  end
  opts.maxit = floor(opts.maxit);

  %% HARDENING CHANGE: required fields
  if ~isfield(data,'bddc') || ~isstruct(data.bddc)
    error('solve_bddc:missing_bddc', 'data.bddc missing (run setup_bddc first).');
  end
  if ~isfield(data.bddc,'R')
    error('solve_bddc:missing_R', 'data.bddc.R missing.');
  end
  if ~isfield(data.bddc,'g_prod')
    error('solve_bddc:missing_g', 'data.bddc.g_prod missing.');
  end

  R = data.bddc.R;

  %% HARDENING CHANGE: validate R
  if ~(isnumeric(R) || issparse(R)) || ndims(R) ~= 2
    error('solve_bddc:bad_R', 'data.bddc.R must be a numeric 2D matrix (dense or sparse).');
  end
  nProd = size(R,1);
  nHat  = size(R,2);

  g_prod = data.bddc.g_prod;
  if ~(isnumeric(g_prod) || issparse(g_prod))
    error('solve_bddc:bad_g', 'data.bddc.g_prod must be numeric.');
  end
  g_prod = g_prod(:); % column

  %% HARDENING CHANGE: dimension checks
  if numel(g_prod) ~= nProd
    error('solve_bddc:dim_mismatch', ...
      'Size mismatch: numel(data.bddc.g_prod)=%d must equal size(R,1)=%d.', numel(g_prod), nProd);
  end

  % Right-hand side in hat space.
  f_hat = R' * g_prod;
  f_hat = f_hat(:); % column (deterministic shape)

  %% HARDENING CHANGE: handle empty hat-space cleanly
  if nHat == 0
    out = struct();
    out.u_hat  = zeros(0,1);
    out.stats  = struct('flag',0,'relres',0,'iter',0,'resvec',[]);
    out.w_prod = R * out.u_hat;
    return;
  end

  %% HARDENING CHANGE: dependency existence checks for clearer errors
  if exist('applyA_hat','file') ~= 2
    error('solve_bddc:missing_dep', 'applyA_hat not found on path.');
  end
  if exist('applyM_bddc','file') ~= 2
    error('solve_bddc:missing_dep', 'applyM_bddc not found on path.');
  end
  if exist('pcg_wrap','file') ~= 2
    error('solve_bddc:missing_dep', 'pcg_wrap not found on path.');
  end

  applyA = @(x) applyA_hat(x, data);
  applyM = @(r) applyM_bddc(r, data);

  %% HARDENING CHANGE: deterministic x0 shape
  x0 = zeros(nHat,1);

  %% HARDENING CHANGE: wrap pcg call to improve failure message context
  try
    [u_hat, stats] = pcg_wrap(applyA, f_hat, opts.tol, opts.maxit, applyM, x0);
  catch err
    error('solve_bddc:pcg_failed', 'pcg_wrap failed in solve_bddc: %s', err.message);
  end

  %% HARDENING CHANGE: normalize returned shape and validate length
  u_hat = u_hat(:);
  if numel(u_hat) ~= nHat
    error('solve_bddc:bad_u_hat', ...
      'pcg_wrap returned u_hat with wrong length: %d (expected %d).', numel(u_hat), nHat);
  end

  out = struct();
  out.u_hat = u_hat;
  out.stats = stats;
  out.w_prod = R * u_hat;
end