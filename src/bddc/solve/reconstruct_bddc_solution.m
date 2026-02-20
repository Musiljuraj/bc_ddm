% ============================================================
% File: src/bddc/solve/reconstruct_bddc_solution.m
% ============================================================
function out = reconstruct_bddc_solution(u_hat, data)
%RECONSTRUCT_BDDC_SOLUTION  Reconstruct global free-DOF solution from u_hat.
%
%   out = reconstruct_bddc_solution(u_hat, data)
%
% Given the hat-space interface solution u_hat, compute:
%   - product interface trace w = R*u_hat
%   - subdomain interior solutions via local back-substitution
%   - averaged global free-DOF solution u_free
%
% Inputs:
%   u_hat : hat-space vector (nHat x 1)
%   data  : struct from build_problem_data(...) extended by setup_bddc(...)
%
% Output fields:
%   out.w_prod : product interface trace (nProd x 1)
%   out.u_free : global free-DOF solution (nFree x 1)
%   out.count  : multiplicity count used in averaging

  if nargin ~= 2
    error('reconstruct_bddc_solution: expected inputs (u_hat, data).');
  end
  if ~isfield(data,'bddc') || ~isfield(data.bddc,'R')
    error('reconstruct_bddc_solution: data.bddc.R missing (run setup_bddc first).');
  end

  R = data.bddc.R;
  w_prod = R * u_hat(:);

  sub = data.sub;
  nSub = numel(sub);
  nFree = data.ddm.nFree;

  u_sum = zeros(nFree, 1);
  cnt = zeros(nFree, 1);

  for i = 1:nSub
    % Interface values for subdomain i (in local gamma ordering).
    pidx = sub(i).prod_idx(:);
    if isempty(pidx)
      continue;
    end
    w_g = w_prod(pidx);

    % Local interior solve: uI = K_II^{-1} (f_I - K_Ig*w_g).
    nI = size(sub(i).K_II, 1);
    if nI > 0
      RII = sub(i).R_II;
      rhsI = sub(i).f_I - sub(i).K_Ig * w_g;
      uI = RII \ (RII' \ rhsI);
    else
      uI = zeros(0,1);
    end

    % Scatter to global free DOFs.
    Iglob = sub(i).glob_I(:);
    Gglob = sub(i).gamma_glob(:);

    if ~isempty(Iglob)
      u_sum(Iglob) = u_sum(Iglob) + uI;
      cnt(Iglob) = cnt(Iglob) + 1;
    end

    if ~isempty(Gglob)
      u_sum(Gglob) = u_sum(Gglob) + w_g;
      cnt(Gglob) = cnt(Gglob) + 1;
    end
  end

  if any(cnt == 0)
    error('reconstruct_bddc_solution: some free DOFs received zero contributions.');
  end

  out = struct();
  out.w_prod = w_prod;
  out.u_free = u_sum ./ cnt;
  out.count = cnt;
end