% ============================================================
% File: src/feti_dp/operators/applyA_lambda.m
% ============================================================
function y = applyA_lambda(x, data)
%APPLYA_LAMBDA  Matrix-free application of the reduced SPD operator A_lambda.
%
% y = A_lambda * x
%
% Operational form (as in your directive):
%   1) t = Bd^T x
%   2) wD = [ \tilde S^{-1} (0, t) ]_Delta
%   3) y = Bd * wD

  x = x(:);
  t = data.BdT * x;                   % packed Delta vector
  out = solve_tildeS([], t, data);    % rc = 0 implicitly
  y = data.Bd * out.wD;
end