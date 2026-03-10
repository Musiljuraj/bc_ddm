% ============================================================
% File: src/bddc/operators/applyM_bddc.m
% ============================================================
function z_hat = applyM_bddc(r_hat, data)
%APPLYM_BDDC Apply the two-level BDDC preconditioner in hat space.
% Thesis link: Chapter 5.4 and Chapter 6.2 (preconditioned BDDC iteration).
% The routine combines scaled distribution, local delta correction, coarse
% correction, and final scaled averaging.
%
% Inputs:
%   r_hat : hat-space residual (nHat x 1)
%   data  : problem struct with fields data.sub, data.primal, data.bddc
%
% Output:
%   z_hat : hat-space vector (nHat x 1)

  if nargin ~= 2
    error('applyM_bddc: expected inputs (r_hat, data).');
  end
  if ~isfield(data,'bddc')
    error('applyM_bddc: data.bddc missing (run setup_bddc first).');
  end

  r_hat = r_hat(:);

  R = data.bddc.R;
  omega = data.bddc.omega;
  Psi = data.bddc.Psi;
  R0 = data.bddc.K0_R;

  % R_D r_hat  (product space)
  r_prod = omega .* (R * r_hat);

  % --- T_sub: local Delta solves ---
  z_sub = zeros(size(r_prod));
  sub = data.sub;
  Sdd_R = data.bddc.Sdd_R;

  for i = 1:numel(sub)
    idx_d = sub(i).prod_idx_d(:);
    if isempty(idx_d)
      continue;
    end
    ri = r_prod(idx_d);
    Rdd = Sdd_R{i};
    zi = Rdd \ (Rdd' \ ri);  % Sdd^{-1} * ri
    z_sub(idx_d) = zi;
  end

  % --- T_0: coarse correction in product space ---
  if isempty(R0)
    z0 = zeros(size(r_prod));
  else
    rhs0 = full(Psi' * r_prod);
    alpha = R0 \ (R0' \ rhs0);
    z0 = Psi * alpha;
  end

  z_prod = z_sub + z0;

  % R_D^T z_prod (hat space)
  z_hat = R' * (omega .* z_prod);
end