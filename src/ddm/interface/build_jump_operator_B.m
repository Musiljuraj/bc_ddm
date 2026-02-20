function [B, meta] = build_jump_operator_B(prod)
%BUILD_JUMP_OPERATOR_B  Construct the jump/constraint operator B enforcing continuity of duplicated interface DOFs.
%
% Link to thesis:
%   Chapter 3.3.3 (jump operator), equations (3.42)–(3.44).
%
% In the product interface space W, each physical interface DOF (assembled index "hat")
% may have k>=2 local copies (product indices). Continuity requires that all copies match.
% We enforce this by constraints:
%   w(p1) - w(pj) = 0,  j=2..k,
% where p1 is a chosen reference copy.
%
% Inputs:
%   prod : product-space bookkeeping from build_product_interface (needs .hat2prod and .nProd/.nHat).
%
% Outputs:
%   B    : sparse constraint matrix of size (m x nProd), with entries in {-1,0,+1}.
%   meta : (optional) struct with fields:
%          .row_hat : (m x 1) assembled interface index associated with each constraint row

  if nargin ~= 1
    error('build_jump_operator_B: expected input (prod).');
  end
  if ~isstruct(prod) || ~isfield(prod,'hat2prod') || ~isfield(prod,'nProd') || ~isfield(prod,'nHat')
    error('build_jump_operator_B: prod must contain hat2prod, nProd, nHat (run build_product_interface).');
  end

  nHat  = prod.nHat;
  nProd = prod.nProd;

  % Count constraints: sum over hat DOFs of (k-1).
  m = 0;
  for h = 1:nHat
    k = numel(prod.hat2prod{h});
    if k > 1
      m = m + (k-1);
    end
  end

  % Triplets for sparse B: each row has exactly two nonzeros (+1 and -1).
  I = zeros(2*m, 1);
  J = zeros(2*m, 1);
  V = zeros(2*m, 1);

  row_hat = zeros(m, 1);

  row = 1;
  ptr = 1;

  %cycle all w-hat elements, take all copies and create constraints for them
  for h = 1:nHat
    copies = prod.hat2prod{h}(:);
    k = numel(copies);
    if k <= 1
      continue;
    end

    p1 = copies(1);
    for j = 2:k  %creating one constraint (star-wise manner, w(p1) with all other w(p's)) pre cycle
      pj = copies(j);

      % Constraint row: w(p1) - w(pj) = 0
      I(ptr)   = row; J(ptr)   = p1; V(ptr)   =  1; ptr = ptr + 1;
      I(ptr)   = row; J(ptr)   = pj; V(ptr)   = -1; ptr = ptr + 1;

      row_hat(row) = h;
      row = row + 1;
    end
  end

  B = sparse(I, J, V, m, nProd);

  if nargout >= 2
    meta = struct();
    meta.row_hat = row_hat;
  end
end
