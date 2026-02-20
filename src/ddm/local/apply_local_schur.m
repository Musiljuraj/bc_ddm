function y = apply_local_schur(si, x)
%APPLY_LOCAL_SCHUR  Matrix-free application of a local Schur complement y = S^(i) x.
%
% Link to thesis:
%   Chapter 3.2.2–3.2.3, equation (3.29).
%
% Given the subdomain block matrices and an interior solver for K_II, compute:
%   y = (K_gg - K_gI * K_II^{-1} * K_Ig) * x
%
% Inputs:
%   si : subdomain struct with fields K_Ig, K_gI, K_gg, R_II (from setup_local_schur).
%   x  : interface vector (length nG for this subdomain).
%
% Output:
%   y  : interface vector y = S^(i) x.
%
% Notes:
% - If the subdomain has no interior DOFs (nI=0), then S^(i) = K_gg.
% - This routine is the core local operator used for S = diag(S_i) on the product interface.

  if nargin ~= 2
    error('apply_local_schur: expected inputs (si, x).');
  end
  if ~isfield(si,'K_gg') || ~isfield(si,'K_Ig') || ~isfield(si,'K_gI')
    error('apply_local_schur: si must contain K_gg, K_Ig, K_gI.');
  end
  if ~isfield(si,'R_II')
    error('apply_local_schur: si must contain R_II (run setup_local_schur first).');
  end

  nG = size(si.K_gg, 1);
  if size(x,1) ~= nG || size(x,2) ~= 1
    error('apply_local_schur: x must be a column vector of length %d.', nG);
  end

  if isempty(si.R_II)
    y = si.K_gg * x;
    return;
  end

  tmp = si.K_Ig * x;                 % nI x 1
  z = si.R_II \ (si.R_II' \ tmp);    % nI x 1
  y   = si.K_gg * x - si.K_gI * z;   % nG x 1
end
