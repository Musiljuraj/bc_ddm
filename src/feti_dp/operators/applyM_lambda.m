% ============================================================
% File: src/feti_dp/operators/applyM_lambda.m
% ============================================================
function z = applyM_lambda(r, data)
%APPLYM_LAMBDA  Baseline FETI-DP preconditioner M_lambda^{-1} (multiplier space).
%
% z = M_lambda^{-1} r
%
% Recipe (directive):
%   t := Bd^T r
%   t := D t
%   u := S_{DeltaDelta} t
%   u := D u
%   z := Bd u
%
% Here S_{DeltaDelta} is the block-diagonal operator with blocks Sdd{i}.

  r = r(:);

  % 1) t = Bd^T r
  t = data.BdT * r;                          % packed Delta

  % 2) t = D t
  t = data.DeltaWeights .* t;

  % 3) u = Sdd * t (blockwise multiply, no solves)
  u = zeros(size(t));
  for i = 1:numel(data.sub)
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    ti = t(rng);
    u(rng) = data.Sdd{i} * ti;
  end

  % 4) u = D u
  u = data.DeltaWeights .* u;

  % 5) z = Bd u
  z = data.Bd * u;
end