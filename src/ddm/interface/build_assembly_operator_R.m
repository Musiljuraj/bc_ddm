function R = build_assembly_operator_R(prod)
%BUILD_ASSEMBLY_OPERATOR_R Build the assembly/distribution operator `R`.
% Thesis link: Chapter 4.3.1 and 4.3.3 (assembled vs. product interface space).
% The operator maps assembled interface vectors to their product-space copies.
%
% Inputs:
%   prod : product-space bookkeeping with fields:
%          .nProd, .nHat, .prod2hat
% Output:
%   R : sparse matrix of size (nProd x nHat)

  if nargin ~= 1
    error('build_assembly_operator_R: expected input (prod).');
  end
  if ~isstruct(prod) || ~isfield(prod,'nProd') || ~isfield(prod,'nHat') || ~isfield(prod,'prod2hat')
    error('build_assembly_operator_R: prod must contain nProd, nHat, prod2hat.');
  end

  nProd = prod.nProd;
  nHat  = prod.nHat;

  % Validate nProd / nHat early for clearer errors
  if ~isnumeric(nProd) || ~isscalar(nProd) || ~isfinite(nProd) || nProd < 0 || nProd ~= fix(nProd)
    error('build_assembly_operator_R: prod.nProd must be a finite nonnegative integer scalar.');
  end
  if ~isnumeric(nHat) || ~isscalar(nHat) || ~isfinite(nHat) || nHat < 0 || nHat ~= fix(nHat)
    error('build_assembly_operator_R: prod.nHat must be a finite nonnegative integer scalar.');
  end

  % Build triplets (I,J,V) for sparse(I,J,V,m,n)
  I = (1:nProd).';         % row indices (1..nProd)
  J = prod.prod2hat(:);    % column index per row

  if ~isnumeric(J)
    error('build_assembly_operator_R: prod.prod2hat must be numeric.');
  end
  if numel(J) ~= nProd
    error('build_assembly_operator_R: prod.prod2hat has wrong length.');
  end

  % Harden: catch NaN/Inf and non-integers explicitly (these may otherwise slip past comparisons)
  if any(~isfinite(J))
    error('build_assembly_operator_R: prod.prod2hat contains NaN/Inf.');
  end
  if any(J ~= fix(J))
    error('build_assembly_operator_R: prod.prod2hat must contain integer indices.');
  end

  if any(J < 1) || any(J > nHat)
    error('build_assembly_operator_R: prod.prod2hat contains out-of-range hat indices.');
  end

  R = sparse(I, J, ones(nProd,1), nProd, nHat);
end