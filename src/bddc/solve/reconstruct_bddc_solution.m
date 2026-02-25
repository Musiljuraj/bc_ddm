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

  %% HARDENING CHANGE: lightweight input validation (types + required top-level fields)
  if ~isnumeric(u_hat) || ~(isvector(u_hat) || isempty(u_hat))
    error('reconstruct_bddc_solution: u_hat must be a numeric vector.');
  end
  if ~isstruct(data)
    error('reconstruct_bddc_solution: data must be a struct.');
  end
  if ~isfield(data,'ddm') || ~isstruct(data.ddm) || ~isfield(data.ddm,'nFree')
    error('reconstruct_bddc_solution: data.ddm.nFree missing.');
  end
  if ~isnumeric(data.ddm.nFree) || ~isscalar(data.ddm.nFree) || ~isfinite(data.ddm.nFree) ...
      || data.ddm.nFree < 0 || abs(data.ddm.nFree - round(data.ddm.nFree)) > 0
    error('reconstruct_bddc_solution: data.ddm.nFree must be a nonnegative integer scalar.');
  end
  if ~isfield(data,'sub') || ~isstruct(data.sub)
    error('reconstruct_bddc_solution: data.sub missing or not a struct array.');
  end

  R = data.bddc.R;

  %% HARDENING CHANGE: explicit size check for R*u_hat (clearer than generic dim mismatch)
  if ~isnumeric(R) || ndims(R) ~= 2
    error('reconstruct_bddc_solution: data.bddc.R must be a numeric 2D matrix.');
  end
  if size(R,2) ~= numel(u_hat)
    error('reconstruct_bddc_solution: size(data.bddc.R,2) must equal numel(u_hat).');
  end

  w_prod = R * u_hat(:);

  sub = data.sub;
  nSub = numel(sub);
  nFree = data.ddm.nFree;

  u_sum = zeros(nFree, 1);
  cnt = zeros(nFree, 1);

  for i = 1:nSub
    %% HARDENING CHANGE: required per-subdomain fields (fail fast with clear message)
    req = {'prod_idx','K_II','R_II','f_I','K_Ig','glob_I','gamma_glob'};
    for k = 1:numel(req)
      if ~isfield(sub(i), req{k})
        error('reconstruct_bddc_solution: sub(%d).%s missing.', i, req{k});
      end
    end

    % Interface values for subdomain i (in local gamma ordering).
    pidx = sub(i).prod_idx(:);

    %% HARDENING CHANGE: validate prod_idx (integer, in range)
    if ~isempty(pidx)
      if ~isnumeric(pidx) || any(~isfinite(pidx)) || any(abs(pidx - round(pidx)) > 0) || any(pidx < 1)
        error('reconstruct_bddc_solution: sub(%d).prod_idx must be a vector of positive integers.', i);
      end
      if any(pidx > numel(w_prod))
        error('reconstruct_bddc_solution: sub(%d).prod_idx out of range for w_prod.', i);
      end
    end

    if isempty(pidx)
      continue;
    end

    w_g = w_prod(pidx);
    nG = numel(w_g);

    % Local interior solve: uI = K_II^{-1} (f_I - K_Ig*w_g).
    KII = sub(i).K_II;
    if ~isnumeric(KII) || ndims(KII) ~= 2
      error('reconstruct_bddc_solution: sub(%d).K_II must be a numeric 2D matrix.', i);
    end
    if size(KII,1) ~= size(KII,2)
      error('reconstruct_bddc_solution: sub(%d).K_II must be square.', i);
    end
    nI = size(KII, 1);

    % Scatter indices (validate early so size checks can refer to them)
    Iglob = sub(i).glob_I(:);
    Gglob = sub(i).gamma_glob(:);

    %% HARDENING CHANGE: gamma_glob must match the local gamma length (if provided)
    if ~isempty(Gglob) && numel(Gglob) ~= nG
      error('reconstruct_bddc_solution: sub(%d).gamma_glob length must match numel(sub(%d).prod_idx).', i, i);
    end

    %% HARDENING CHANGE: validate global index vectors (integer, in range)
    if ~isempty(Iglob)
      if ~isnumeric(Iglob) || any(~isfinite(Iglob)) || any(abs(Iglob - round(Iglob)) > 0) || any(Iglob < 1) || any(Iglob > nFree)
        error('reconstruct_bddc_solution: sub(%d).glob_I must be integer indices in [1..nFree].', i);
      end
    end
    if ~isempty(Gglob)
      if ~isnumeric(Gglob) || any(~isfinite(Gglob)) || any(abs(Gglob - round(Gglob)) > 0) || any(Gglob < 1) || any(Gglob > nFree)
        error('reconstruct_bddc_solution: sub(%d).gamma_glob must be integer indices in [1..nFree].', i);
      end
    end

    if nI > 0
      RII = sub(i).R_II;
      fI  = sub(i).f_I;
      KIg = sub(i).K_Ig;

      %% HARDENING CHANGE: local dimension consistency checks (cheap)
      if ~isnumeric(RII) || ~isequal(size(RII), [nI nI])
        error('reconstruct_bddc_solution: sub(%d).R_II must be %dx%d.', i, nI, nI);
      end
      if ~isnumeric(fI) || numel(fI) ~= nI
        error('reconstruct_bddc_solution: sub(%d).f_I must have length %d.', i, nI);
      end
      if ~isnumeric(KIg) || ndims(KIg) ~= 2 || size(KIg,1) ~= nI || size(KIg,2) ~= nG
        error('reconstruct_bddc_solution: sub(%d).K_Ig must be %dx%d.', i, nI, nG);
      end
      if ~isempty(Iglob) && numel(Iglob) ~= nI
        error('reconstruct_bddc_solution: sub(%d).glob_I length must equal nI=%d.', i, nI);
      end

      rhsI = fI(:) - KIg * w_g;

      %% HARDENING CHANGE: catch solve failures with subdomain context
      try
        uI = RII \ (RII' \ rhsI);
      catch ME
        error('reconstruct_bddc_solution: interior solve failed on subdomain %d (%s).', i, ME.message);
      end
    else
      uI = zeros(0,1);
    end

    %% HARDENING CHANGE: ensure scatter sizes match computed vectors
    if ~isempty(Iglob) && numel(uI) ~= numel(Iglob)
      error('reconstruct_bddc_solution: sub(%d) interior vector size mismatch (uI vs glob_I).', i);
    end
    if ~isempty(Gglob) && numel(w_g) ~= numel(Gglob)
      error('reconstruct_bddc_solution: sub(%d) gamma vector size mismatch (w_g vs gamma_glob).', i);
    end

    % Scatter to global free DOFs.
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
    %% HARDENING CHANGE: slightly clearer diagnostic
    z = find(cnt == 0);
    error('reconstruct_bddc_solution: some free DOFs received zero contributions (e.g. first is %d).', z(1));
  end

  out = struct();
  out.w_prod = w_prod;
  out.u_free = u_sum ./ cnt;
  out.count = cnt;
end