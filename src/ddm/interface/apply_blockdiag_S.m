function y = apply_blockdiag_S(sub, w)
%APPLY_BLOCKDIAG_S  Apply the block-diagonal interface operator S = diag(S_i) in the product space.
%
% Link to thesis:
%   Chapter 3.4.1, equation (3.54).
%
% After interior elimination, each subdomain contributes a local Schur complement S^(i)
% acting on its local interface DOFs. In the product (duplicated) space, the global
% operator is block-diagonal:
%     S = diag(S^(1), ..., S^(N)).
%
% Inputs:
%   sub : subdomain array with fields:
%         .prod_idx (product indices for this subdomain's interface DOFs)
%         and the local Schur data required by apply_local_schur (K blocks + R_II).
%   w   : product interface vector (length nProd).
%
% Output:
%   y   : product interface vector y = S*w.

  if nargin ~= 2
    error('apply_blockdiag_S: expected inputs (sub, w).');
  end
  if ~isvector(w) || size(w,2) ~= 1
    error('apply_blockdiag_S: w must be a column vector.');
  end

  nProd = numel(w);
  y = zeros(nProd, 1);

  nSub = numel(sub);
  for i = 1:nSub
    if ~isfield(sub(i),'prod_idx')
      error('apply_blockdiag_S: sub(%d) missing prod_idx (run build_product_interface).', i);
    end
    idx = sub(i).prod_idx(:);
    if isempty(idx)
      continue;
    end
    if any(idx < 1) || any(idx > nProd)
      error('apply_blockdiag_S: sub(%d).prod_idx contains out-of-range indices.', i);
    end

    wi = w(idx);
    yi = apply_local_schur(sub(i), wi);
    y(idx) = yi;
  end
end
