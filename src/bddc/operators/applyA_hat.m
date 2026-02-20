% ============================================================
% File: src/bddc/operators/applyA_hat.m
% ============================================================
function y = applyA_hat(x, data)
%APPLYA_HAT  Apply hat-space operator A = R^T S R.
%
%   y = applyA_hat(x, data)
%
% Inputs:
%   x    : hat-space vector (nHat x 1)
%   data : problem struct with fields:
%          - sub (from build_problem_data)
%          - bddc.R (product assembly operator)
%
% Output:
%   y : hat-space vector (nHat x 1)

  if nargin ~= 2
    error('applyA_hat: expected inputs (x, data).');
  end
  if ~isfield(data,'bddc') || ~isfield(data.bddc,'R')
    error('applyA_hat: data.bddc.R missing (run setup_bddc first).');
  end

  x = x(:);
  R = data.bddc.R;

  % product lift
  w = R * x;

  % apply block-diagonal Schur in product space
  Sw = apply_blockdiag_S(data.sub, w);

  % assemble back to hat space
  y = R' * Sw;
end