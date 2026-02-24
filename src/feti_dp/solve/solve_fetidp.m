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

  narginchk(1,3); % ADDED

  if nargin < 2 || isempty(tol);   tol = 1e-8; end
  if nargin < 3 || isempty(maxit); maxit = 500; end

  % -----------------------------
  % Basic input validation
  % -----------------------------
  if ~isstruct(data) % ADDED
    error('solve_fetidp:invalidData', 'data must be a struct.'); % ADDED
  end

  % Required top-level fields (minimal contract) % ADDED
  req = {'sub','delta_range','nDeltaProd','primal','Bd','nLambda'}; % ADDED
  for k = 1:numel(req) % ADDED
    if ~isfield(data, req{k}) % ADDED
      error('solve_fetidp:missingField', 'data.%s is required.', req{k}); % ADDED
    end
  end

  if ~isstruct(data.primal) || ~isfield(data.primal,'nC') % ADDED
    error('solve_fetidp:invalidData', 'data.primal.nC is required.'); % ADDED
  end

  if ~(isnumeric(tol) && isscalar(tol) && isreal(tol) && isfinite(tol) && tol > 0) % ADDED
    error('solve_fetidp:invalidTol', 'tol must be a real, finite, positive scalar.'); % ADDED
  end
  if ~(isnumeric(maxit) && isscalar(maxit) && isreal(maxit) && isfinite(maxit) && maxit >= 1 && maxit == floor(maxit)) % ADDED
    error('solve_fetidp:invalidMaxit', 'maxit must be a real, finite, positive integer scalar.'); % ADDED
  end

  if ~(isnumeric(data.nLambda) && isscalar(data.nLambda) && isfinite(data.nLambda) && data.nLambda >= 1 && data.nLambda == floor(data.nLambda)) % ADDED
    error('solve_fetidp:invalidData', 'data.nLambda must be a positive integer scalar.'); % ADDED
  end
  if ~(isnumeric(data.nDeltaProd) && isscalar(data.nDeltaProd) && isfinite(data.nDeltaProd) && data.nDeltaProd >= 0 && data.nDeltaProd == floor(data.nDeltaProd)) % ADDED
    error('solve_fetidp:invalidData', 'data.nDeltaProd must be a nonnegative integer scalar.'); % ADDED
  end
  if ~(isnumeric(data.primal.nC) && isscalar(data.primal.nC) && isfinite(data.primal.nC) && data.primal.nC >= 0 && data.primal.nC == floor(data.primal.nC)) % ADDED
    error('solve_fetidp:invalidData', 'data.primal.nC must be a nonnegative integer scalar.'); % ADDED
  end

  if ~(isnumeric(data.Bd) && isreal(data.Bd) && all(isfinite(data.Bd(:)))) % ADDED
    error('solve_fetidp:invalidData', 'data.Bd must be a real, finite numeric matrix.'); % ADDED
  end
  if size(data.Bd,1) ~= data.nLambda % ADDED
    error('solve_fetidp:dimMismatch', 'size(data.Bd,1) must equal data.nLambda.'); % ADDED
  end
  if size(data.Bd,2) ~= data.nDeltaProd % ADDED
    error('solve_fetidp:dimMismatch', 'size(data.Bd,2) must equal data.nDeltaProd.'); % ADDED
  end

  if ~isstruct(data.sub) % ADDED
    error('solve_fetidp:invalidData', 'data.sub must be a struct array.'); % ADDED
  end
  sub  = data.sub;
  nSub = numel(sub);

  if ~(iscell(data.delta_range) && numel(data.delta_range) == nSub) % ADDED
    error('solve_fetidp:invalidData', 'data.delta_range must be a cell array with one entry per subdomain.'); % ADDED
  end

  % Matrix-free operator and preconditioner in multiplier space
  applyA = @(x) applyA_lambda(x, data);
  applyM = @(r) applyM_lambda(r, data);

  % ------------------------------------------------------------
  % Build RHS pieces from local condensed interface RHS sub(i).g
  %
  % gD : packed Delta RHS (length nDeltaProd)
  % gC : assembled coarse/primal RHS in global coarse indexing (length nC)
  % ------------------------------------------------------------
  gD = zeros(data.nDeltaProd, 1);
  nC = data.primal.nC;
  gC = zeros(nC, 1);

  for i = 1:nSub
    % Validate required subfields (minimal) % ADDED
    reqs = {'g','idx_d','idx_c','c_ids'}; % ADDED
    for k = 1:numel(reqs) % ADDED
      if ~isfield(sub(i), reqs{k}) % ADDED
        error('solve_fetidp:missingField', 'data.sub(%d).%s is required.', i, reqs{k}); % ADDED
      end
    end

    gi = sub(i).g(:);
    if ~(isnumeric(gi) && isreal(gi) && all(isfinite(gi))) % ADDED
      error('solve_fetidp:invalidData', 'data.sub(%d).g must be a real, finite numeric vector.', i); % ADDED
    end

    idx_d = sub(i).idx_d(:);
    idx_c = sub(i).idx_c(:);
    c_ids = sub(i).c_ids(:);

    if ~(isnumeric(idx_d) && isreal(idx_d) && all(isfinite(idx_d)) && all(idx_d == floor(idx_d))) % ADDED
      error('solve_fetidp:invalidIndex', 'data.sub(%d).idx_d must be integer-valued.', i); % ADDED
    end
    if ~(isnumeric(idx_c) && isreal(idx_c) && all(isfinite(idx_c)) && all(idx_c == floor(idx_c))) % ADDED
      error('solve_fetidp:invalidIndex', 'data.sub(%d).idx_c must be integer-valued.', i); % ADDED
    end
    if ~(isnumeric(c_ids) && isreal(c_ids) && all(isfinite(c_ids)) && all(c_ids == floor(c_ids))) % ADDED
      error('solve_fetidp:invalidIndex', 'data.sub(%d).c_ids must be integer-valued.', i); % ADDED
    end

    if ~isempty(idx_d) % ADDED
      if any(idx_d < 1) || any(idx_d > numel(gi)) % ADDED
        error('solve_fetidp:indexRange', 'data.sub(%d).idx_d out of range for sub(%d).g.', i, i); % ADDED
      end
    end
    if ~isempty(idx_c) % ADDED
      if any(idx_c < 1) || any(idx_c > numel(gi)) % ADDED
        error('solve_fetidp:indexRange', 'data.sub(%d).idx_c out of range for sub(%d).g.', i, i); % ADDED
      end
    end
    if ~isempty(c_ids) % ADDED
      if any(c_ids < 1) || any(c_ids > nC) % ADDED
        error('solve_fetidp:indexRange', 'data.sub(%d).c_ids out of range 1..nC.', i); % ADDED
      end
      if numel(c_ids) ~= numel(idx_c) % ADDED
        error('solve_fetidp:dimMismatch', 'sub(%d): numel(c_ids) must match numel(idx_c).', i); % ADDED
      end
    end

    % Delta part (packed using data.delta_range{i})
    rng = data.delta_range{i};
    if ~isempty(rng)
      rng = rng(:);
      if ~(isnumeric(rng) && isreal(rng) && all(isfinite(rng)) && all(rng == floor(rng))) % ADDED
        error('solve_fetidp:invalidIndex', 'data.delta_range{%d} must be integer-valued.', i); % ADDED
      end
      if any(rng < 1) || any(rng > data.nDeltaProd) % ADDED
        error('solve_fetidp:indexRange', 'data.delta_range{%d} out of range 1..nDeltaProd.', i); % ADDED
      end
      if numel(rng) ~= numel(idx_d) % ADDED
        error('solve_fetidp:dimMismatch', 'delta_range{%d} length must match numel(sub(%d).idx_d).', i, i); % ADDED
      end
      gD(rng) = gi(idx_d);
    end

    % Primal/coarse part (assembled using sub(i).c_ids)
    if ~isempty(c_ids)
      gc_i = gi(idx_c);
      gC(c_ids) = gC(c_ids) + gc_i;
    end
  end

  % ------------------------------------------------------------
  % Form multiplier-space RHS:
  %   b = Bd * wD0,  where (u_c0, wD0) = \tilde S^{-1} (gC, gD)
  % ------------------------------------------------------------
  out0 = solve_tildeS(gC, gD, data);
  if ~(isstruct(out0) && isfield(out0,'wD')) % ADDED
    error('solve_fetidp:invalidOutput', 'solve_tildeS must return a struct with field wD.'); % ADDED
  end
  wD0 = out0.wD(:); % ADDED
  if ~(isnumeric(wD0) && isreal(wD0) && all(isfinite(wD0))) % ADDED
    error('solve_fetidp:invalidOutput', 'solve_tildeS output wD must be real and finite.'); % ADDED
  end
  if numel(wD0) ~= data.nDeltaProd % ADDED
    error('solve_fetidp:dimMismatch', 'solve_tildeS output wD length must equal data.nDeltaProd.'); % ADDED
  end

  b = data.Bd * wD0;
  if ~(isvector(b) && numel(b) == data.nLambda) % ADDED
    error('solve_fetidp:dimMismatch', 'Multiplier RHS b must have length data.nLambda.'); % ADDED
  end

  % Initial guess
  x0 = zeros(data.nLambda, 1);

  % PCG
  [lambda, stats] = pcg_wrap(applyA, b, tol, maxit, applyM, x0);

end