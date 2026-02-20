function omega = multiplicity_scaling(prod)
%MULTIPLICITY_SCALING  Compute multiplicity scaling weights on product interface DOFs.
%
%   omega = multiplicity_scaling(prod)
%
% For each assembled interface DOF (hat-index), there are k product copies.
% This returns product-space weights omega such that each copy gets weight 1/k.
%
% Required inputs:
%   prod.nProd
%   prod.hat2prod  (cell array; hat2prod{h} lists product indices)
%
% Output:
%   omega (nProd x 1), omega(p) = 1/k in each hat group of size k.

  if nargin ~= 1
    error('multiplicity_scaling: expected exactly one input (prod).');
  end
  if ~isstruct(prod) || ~isfield(prod,'nProd') || ~isfield(prod,'hat2prod')
    error('multiplicity_scaling: prod must contain fields nProd and hat2prod.');
  end

  nProd = prod.nProd;
  hat2prod = prod.hat2prod;

  if ~(isscalar(nProd) && nProd == floor(nProd) && nProd >= 0)
    error('multiplicity_scaling: prod.nProd must be an integer >= 0.');
  end
  if ~iscell(hat2prod)
    error('multiplicity_scaling: prod.hat2prod must be a cell array.');
  end

  omega = zeros(nProd, 1);

  for h = 1:numel(hat2prod)
    idx = hat2prod{h};
    if isempty(idx)
      continue;
    end
    idx = idx(:);
    k = numel(idx);
    omega(idx) = 1.0 / k;
  end

  if any(omega == 0)
    error('multiplicity_scaling: some product indices received zero weight (check prod.hat2prod).');
  end
end