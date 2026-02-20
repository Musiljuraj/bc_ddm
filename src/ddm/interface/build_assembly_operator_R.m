function R = build_assembly_operator_R(prod)
%BUILD_ASSEMBLY_OPERATOR_R  Build the assembly (distribution) operator R: \hat{W} -> W.
%
% Link to thesis:
%   Chapter 3.3.3, equation (3.48).
%
% For an assembled interface vector \hat{w} (single-valued per physical interface DOF),
% the product representation w = R * \hat{w} duplicates values to all local copies.
% R has entries 0/1 and exactly one '1' per product row:
%   w(p) = \hat{w}( prod.prod2hat(p) ).
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

  %The code builds triplets (I,J,V) and calls Octave’s sparse(I,J,V,m,n) constructor:
  I = (1:nProd).'; %create row indices for each product DOF row (1..nProd).
  J = prod.prod2hat(:); %column index per row, telling which hat DOF each product DOF copies from

  if numel(J) ~= nProd
    error('build_assembly_operator_R: prod.prod2hat has wrong length.');
  end
  if any(J < 1) || any(J > nHat)
    error('build_assembly_operator_R: prod.prod2hat contains out-of-range hat indices.');
  end

  R = sparse(I, J, ones(nProd,1), nProd, nHat); %build the 0/1 distribution matrix with one “1” per product DOF row.
end
