function y = apply_local_schur(si, x)
%APPLY_LOCAL_SCHUR Apply one local Schur complement in matrix-free form.
% Thesis link: Chapter 4.2.3 (local Schur complement operator).
% For a given subdomain block structure, the routine computes `y = S^(i) x`
% without forming `S^(i)` explicitly.

  if nargin ~= 2
    error('apply_local_schur:BadNargin', 'Expected inputs (si, x).');
  end
  if ~isstruct(si)
    error('apply_local_schur:BadInput', 'si must be a struct.');
  end
  req = {'K_gg','K_Ig','K_gI','R_II'};
  for k = 1:numel(req)
    if ~isfield(si, req{k})
      error('apply_local_schur:MissingField', 'si must contain field "%s".', req{k});
    end
  end
  if ~isnumeric(x)
    error('apply_local_schur:BadInput', 'x must be numeric.');
  end

  Kgg = si.K_gg; KIg = si.K_Ig; KgI = si.K_gI; R = si.R_II;

  % Basic block dimension checks
  nG = size(Kgg,1);
  if size(Kgg,2) ~= nG
    error('apply_local_schur:DimMismatch', 'K_gg must be square (nG x nG).');
  end
  if size(KIg,2) ~= nG
    error('apply_local_schur:DimMismatch', 'K_Ig must have size (nI x nG).');
  end
  nI = size(KIg,1);
  if size(KgI,1) ~= nG || size(KgI,2) ~= nI
    error('apply_local_schur:DimMismatch', 'K_gI must have size (nG x nI).');
  end

  % x shape: allow multiple RHS
  if size(x,1) ~= nG
    error('apply_local_schur:DimMismatch', 'x must have %d rows.', nG);
  end

  if isempty(R)
    if nI ~= 0
      error('apply_local_schur:InconsistentState', ...
            'R_II is empty but nI=%d (expected nI=0).', nI);
    end
    y = Kgg * x;
    return;
  end

  if size(R,1) ~= nI || size(R,2) ~= nI
    error('apply_local_schur:DimMismatch', 'R_II must be size (nI x nI).');
  end

  % Apply Schur action
  tmp = KIg * x;

  % Robust to whether R is upper or lower triangular
  if istriu(R)
    z = R \ (R' \ tmp);      % for K_II = R'*R
  elseif istril(R)
    z = R' \ (R \ tmp);      % for K_II = R*R'
  else
    error('apply_local_schur:BadFactor', 'R_II must be triangular (upper or lower Cholesky factor).');
  end

  y = Kgg * x - KgI * z;
end