% ============================================================
% File: src/feti_dp/operators/solve_tildeS.m
% ============================================================
function out = solve_tildeS(rc, rD, data)
%SOLVE_TILDES  Apply \widetilde{S}^{-1} to RHS split into primal and delta parts.
%
% This implements the core reduced SPD solve in product interface space:
% given (rc, rD), returns w_c and w_D such that w = \tilde S^{-1} [rc; rD].
%
% Packed representations:
% - rD is a packed Delta vector (length nDeltaProd) consistent with primal.prod2delta ordering.
% - rc is given implicitly in assembled coarse coordinates (length nC) OR as zeros.
%   In this implementation we treat rc as "product-primal RHS assembled into coarse ids"
%   because only A_i^T(...) appears in the coarse rhs formula.
%
% Output "out" struct:
%   out.u_c        : coarse vector (nC)
%   out.wc_local{i}: local primal values on subdomain i (length nC_i)
%   out.wD         : packed Delta vector (nDeltaProd)

  sub = data.sub;
  nSub = numel(sub);

  nC = data.primal.nC;
  if nargin < 1 || isempty(rc); rc = zeros(nC,1); end
  if nargin < 2 || isempty(rD); rD = zeros(data.nDeltaProd,1); end

  rc = rc(:);
  rD = rD(:);

  % ------------------------------------------------------------
  % Step 1: yD = Sdd^{-1} rD (blockwise)
  % ------------------------------------------------------------
  yD = zeros(size(rD));
  yD_local = cell(nSub,1);

  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      yD_local{i} = zeros(0,1);
      continue;
    end
    ri = rD(rng);
    R = data.Sdd_R{i};
    yi = R \ (R' \ ri);
    yD(rng) = yi;
    yD_local{i} = yi;
  end

  % ------------------------------------------------------------
  % Step 2: coarse rhs: rhs_c = rc - sum_i Scd_i * yD_i  assembled into coarse ids
  % ------------------------------------------------------------
  rhs_c = rc;
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      continue;
    end
    yi = yD_local{i};
    if isempty(yi)
      continue;
    end
    rhs_c(c_ids) = rhs_c(c_ids) - data.Scd{i} * yi;
  end

  % ------------------------------------------------------------
  % Step 3: solve u_c = Kcc^{-1} rhs_c
  % ------------------------------------------------------------
  if nC > 0
    Rcc = data.Kcc_R;
    u_c = Rcc \ (Rcc' \ rhs_c);
  else
    u_c = zeros(0,1);
  end

  % ------------------------------------------------------------
  % Step 4: w_c local on each subdomain (just u_c restricted by c_ids)
  % ------------------------------------------------------------
  wc_local = cell(nSub,1);
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      wc_local{i} = zeros(0,1);
    else
      wc_local{i} = u_c(c_ids);
    end
  end

  % ------------------------------------------------------------
  % Step 5: wD = yD - Sdd^{-1}( Sdc_i * wci ) blockwise
  % ------------------------------------------------------------
  wD = yD;

  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    wci = wc_local{i};
    if isempty(wci)
      % no primal on this subdomain -> correction is zero
      continue;
    end
    corr = data.Sdc{i} * wci;
    R = data.Sdd_R{i};
    z = R \ (R' \ corr);
    wD(rng) = wD(rng) - z;
  end

  out = struct();
  out.u_c = u_c;
  out.wc_local = wc_local;
  out.wD = wD;
end