% ============================================================
% File: src/feti_dp/solve/solve_fetidp.m
% ============================================================
function [lambda, stats] = solve_fetidp(data, tol, maxit)
%SOLVE_FETIDP  Run PCG in multiplier space on the reduced SPD operator A_lambda.
%
% Solves:  A_lambda * lambda = b
% where:   A_lambda = Bd * [ \tilde S^{-1} (0, Bd^T * . ) ]_Delta
%
% IMPORTANT CHANGE (vs previous version):
%   The multiplier RHS b must be formed using the FULL condensed interface RHS,
%   i.e. BOTH primal/coarse part gC and delta part gD:
%
%     b = Bd * [ \tilde S^{-1} (gC, gD) ]_Delta
%
%   Previously the code used rc=0:
%     b = Bd * [ \tilde S^{-1} (0, gD) ]_Delta
%   This is generally incorrect when nC>0 and gC is nonzero, because PCG then
%   solves a different dual problem than the one used in reconstruction, leading
%   to a noticeable constraint residual ||Bd*wD|| after reconstruction.
%
% Inputs:
%   data  : struct produced by build_problem_data + setup_fetidp
%   tol   : PCG tolerance (default 1e-8)
%   maxit : PCG max iterations (default 500)
%
% Outputs:
%   lambda : reduced multiplier vector (length data.nLambda)
%   stats  : pcg_wrap stats (flag, relres, iter, resvec)

  if nargin < 2 || isempty(tol);   tol = 1e-8; end
  if nargin < 3 || isempty(maxit); maxit = 500; end

  % Matrix-free operator and preconditioner in multiplier space
  applyA = @(x) applyA_lambda(x, data);
  applyM = @(r) applyM_lambda(r, data);

  sub  = data.sub;
  nSub = numel(sub);

  % ------------------------------------------------------------
  % Build RHS pieces from local condensed interface RHS sub(i).g
  %
  % gD : packed Delta RHS (length nDeltaProd)
  % gC : assembled coarse/primal RHS in global coarse indexing (length nC)
  %
  % NOTE: gC is NOT the same as data.hc.
  % - data.hc is the already assembled "coarse condensed RHS"
  %     hc = sum_i ( gc_i - Scd_i * Sdd_i^{-1} * gd_i )
  % - Here we need the raw assembled primal component gC = sum_i gc_i,
  %   because solve_tildeS internally performs the "-Scd*Sdd^{-1}*gD" coupling
  %   when forming the coarse rhs.
  % ------------------------------------------------------------
  gD = zeros(data.nDeltaProd, 1);

  nC = data.primal.nC;
  gC = zeros(nC, 1);

  for i = 1:nSub
    gi = sub(i).g(:);

    % Delta part (packed using data.delta_range{i})
    rng = data.delta_range{i};
    if ~isempty(rng)
      gD(rng) = gi(sub(i).idx_d(:));
    end

    % Primal/coarse part (assembled using sub(i).c_ids)
    c_ids = sub(i).c_ids(:);
    if ~isempty(c_ids)
      gc_i = gi(sub(i).idx_c(:));
      % Assemble into global coarse coordinates
      gC(c_ids) = gC(c_ids) + gc_i;
    end
  end

  % ------------------------------------------------------------
  % Form multiplier-space RHS:
  %   b = Bd * wD0,  where (u_c0, wD0) = \tilde S^{-1} (gC, gD)
  % ------------------------------------------------------------
  out0 = solve_tildeS(gC, gD, data);
  b = data.Bd * out0.wD;

  % Initial guess
  x0 = zeros(data.nLambda, 1);

  % PCG
  [lambda, stats] = pcg_wrap(applyA, b, tol, maxit, applyM, x0);

end