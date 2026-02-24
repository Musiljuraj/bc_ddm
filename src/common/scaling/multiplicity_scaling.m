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

  % CHANGED: stronger validation of nProd (numeric, real, finite, integer >= 0; reject logical)
  if ~(isnumeric(nProd) && ~islogical(nProd) && isreal(nProd) && isfinite(nProd) && ...
       isscalar(nProd) && nProd == floor(nProd) && nProd >= 0)
    error('multiplicity_scaling: prod.nProd must be a real, finite integer >= 0.');
  end

  if ~iscell(hat2prod)
    error('multiplicity_scaling: prod.hat2prod must be a cell array.');
  end

  omega = zeros(nProd, 1);

  % ADDED: track coverage to detect overlaps and to make the final check explicit
  seen = false(nProd, 1);  % seen(p)=true iff product index p assigned already

  for h = 1:numel(hat2prod)
    idx = hat2prod{h};

    if isempty(idx)
      continue;
    end

    % ADDED: validate hat2prod{h} element type/shape
    if ~(isnumeric(idx) && ~islogical(idx) && isreal(idx))
      error('multiplicity_scaling: hat2prod{%d} must be a real, non-logical numeric vector of indices.', h);
    end
    if ~isvector(idx)
      error('multiplicity_scaling: hat2prod{%d} must be a vector of indices.', h);
    end
    if any(~isfinite(idx(:)))
      error('multiplicity_scaling: hat2prod{%d} contains NaN/Inf.', h);
    end

    idx = idx(:);

    % ADDED: integer-valued index requirement
    if any(idx ~= floor(idx))
      error('multiplicity_scaling: hat2prod{%d} must contain integer-valued indices.', h);
    end

    % ADDED: range check (1..nProd)
    if any(idx < 1) || any(idx > nProd)
      error('multiplicity_scaling: hat2prod{%d} contains indices outside 1..nProd.', h);
    end

    % ADDED: reject duplicates within a hat group
    if numel(unique(idx)) ~= numel(idx)
      error('multiplicity_scaling: hat2prod{%d} contains duplicate indices.', h);
    end

    % FIXED: reject overlaps across hat groups (previously silent overwrite)
    if any(seen(idx))
      error('multiplicity_scaling: a product index appears in multiple hat groups (overlap).');
    end

    k = numel(idx);
    omega(idx) = 1.0 / k;
    seen(idx) = true;  % ADDED
  end

  % CHANGED: coverage check expressed via seen mask (equivalent to omega==0 under our new rules)
  if any(~seen)
    error('multiplicity_scaling: some product indices received zero weight (check prod.hat2prod).');
  end
end