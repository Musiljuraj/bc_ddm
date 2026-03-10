function [x, stats] = pcg_wrap(applyA, b, tol, maxit, applyM, x0)
%PCG_WRAP Thin wrapper around Octave's `pcg` for matrix-free SPD solves.
% Thesis link: Chapter 5.2 (common solver infrastructure).
% The routine standardizes callback-based PCG calls and returns a compact
% diagnostics structure used in later experiments.
%
% Inputs:
%   applyA : function handle, y = applyA(x). Must represent an SPD operator.
%   b      : right-hand side vector.
%   tol    : relative residual tolerance passed to pcg.
%   maxit  : maximum number of iterations passed to pcg.
%   applyM : optional preconditioner handle z = applyM(r), approximating A\r.
%            If [] or omitted, no preconditioner is used.
%   x0     : optional initial guess. If [] or omitted, zeros(size(b)) is used.
%
% Outputs:
%   x     : approximate solution.
%   stats : struct with fields:
%           - flag   : pcg exit flag
%           - relres : achieved relative residual
%           - iter   : iteration count
%           - resvec : residual history (length iter+1)


  if nargin < 4
    error('pcg_wrap: expected at least 4 inputs (applyA, b, tol, maxit).');
  end

  if ~isa(applyA, 'function_handle')
    error('pcg_wrap: applyA must be a function handle.');
  end

  if ~isvector(b)
    error('pcg_wrap: b must be a vector.');
  end
  b = b(:);

  if ~(isscalar(tol) && isfinite(tol) && tol > 0)
    error('pcg_wrap: tol must be a positive finite scalar.');
  end

  if ~(isscalar(maxit) && isfinite(maxit) && maxit == floor(maxit) && maxit >= 0)
    error('pcg_wrap: maxit must be an integer >= 0.');
  end

  if nargin < 5
    applyM = [];
  end

  if nargin < 6
    x0 = [];
  end

  if isempty(x0)
    x0 = zeros(size(b));
  else
    if ~isvector(x0)
      error('pcg_wrap: x0 must be a vector (or []).');
    end
    x0 = x0(:);
    if numel(x0) ~= numel(b)
      error('pcg_wrap: x0 must have the same length as b.');
    end
  end

  % Call Octave built-in pcg. Use function-handle form.
  % Preconditioner:
  %   - if applyM is empty -> pass []
  %   - else -> pass applyM as M1 (factor of the preconditioner)
  if isempty(applyM)
    [x, flag, relres, iter, resvec] = pcg(applyA, b, tol, maxit, [], [], x0);
  else
    if ~isa(applyM, 'function_handle')
      error('pcg_wrap: applyM must be a function handle or [].');
    end
    [x, flag, relres, iter, resvec] = pcg(applyA, b, tol, maxit, applyM, [], x0);
  end

  stats = struct();
  stats.flag   = flag;
  stats.relres = relres;
  stats.iter   = iter;
  stats.resvec = resvec;
end