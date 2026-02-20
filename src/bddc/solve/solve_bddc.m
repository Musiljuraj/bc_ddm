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

  if nargin < 2
    opts = struct();
  end
  if ~isfield(opts,'tol');   opts.tol = 1e-8; end
  if ~isfield(opts,'maxit'); opts.maxit = 500; end

  if ~isfield(data,'bddc')
    error('solve_bddc: data.bddc missing (run setup_bddc first).');
  end

  R = data.bddc.R;

  % Right-hand side in hat space.
  f_hat = R' * data.bddc.g_prod;

  applyA = @(x) applyA_hat(x, data);
  applyM = @(r) applyM_bddc(r, data);

  x0 = zeros(size(f_hat));
  [u_hat, stats] = pcg_wrap(applyA, f_hat, opts.tol, opts.maxit, applyM, x0);

  out = struct();
  out.u_hat = u_hat;
  out.stats = stats;
  out.w_prod = R * u_hat;
end