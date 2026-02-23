function [sub, prod] = build_product_interface(sub, iface)
%BUILD_PRODUCT_INTERFACE  Construct product (duplicated) interface indexing and bookkeeping.
%
% Link to thesis:
%   Chapter 3.3.2 (product interface space), equations (3.45)–(3.47).
%
% This routine constructs the product interface space W = Π_i W_i by assigning
% a unique "product" index to every local interface DOF on every subdomain.
% Physical interface DOFs (assembled interface space \hat{W}) are represented
% by iface.glob2hat.
%
% Inputs:
%   sub   : subdomain array with field .glob_G (global free DOFs on interface)
%   iface : interface bookkeeping (needs .glob2hat and .nHat)
%
% Outputs:
%   sub  : updated with fields:
%     .gamma_glob : global free DOFs on interface (column)
%     .gamma_hat  : assembled interface indices (column, 1..nHat)
%     .prod_idx   : product indices for this subdomain's interface vector (contiguous segment)
%
%   prod : product-space bookkeeping:
%     .nHat      : number of assembled interface DOFs
%     .nProd     : dimension of product space (sum_i |glob_G^(i)|)
%     .prod2sub  : (nProd x 1) subdomain id for each product DOF
%     .prod2hat  : (nProd x 1) assembled interface index for each product DOF
%     .hat2prod  : cell (nHat x 1): hat2prod{k} lists product copies of hat DOF k

  % CHANGED: keep strict arg-count check (also supports unit test)
  if nargin ~= 2
    error('build_product_interface: expected inputs (sub, iface).');
  end

  % CHANGED: validate sub container
  if ~isstruct(sub)
    error('build_product_interface: sub must be a struct array.');
  end

  % CHANGED: validate iface required fields early
  if ~isstruct(iface) || ~isfield(iface,'glob2hat') || ~isfield(iface,'nHat')
    error('build_product_interface: iface must contain glob2hat and nHat (run identify_interface_dofs).');
  end

  % CHANGED: validate iface.nHat (numeric, real, finite, integer scalar, >=0)
  nHat = iface.nHat;
  if ~isnumeric(nHat) || ~isreal(nHat) || ~isfinite(nHat) || ~isscalar(nHat) || (nHat ~= fix(nHat)) || (nHat < 0)
    error('build_product_interface: iface.nHat must be a real, finite, integer scalar >= 0.');
  end
  nHat = double(nHat); % normalize type for downstream

  % CHANGED: validate iface.glob2hat (numeric, real, finite, vector; reject char/logical)
  g2h = iface.glob2hat;
  if ischar(g2h) || islogical(g2h) || (exist('isstring','file')==2 && isstring(g2h))
    error('build_product_interface: iface.glob2hat must be numeric (char/logical/string not allowed).');
  end
  if ~isnumeric(g2h) || ~isreal(g2h) || ~isvector(g2h) || (~isempty(g2h) && any(~isfinite(g2h(:))))
    error('build_product_interface: iface.glob2hat must be a real, finite numeric vector (possibly empty).');
  end
  if ~isempty(g2h)
    if any(g2h(:) ~= fix(g2h(:)))
      error('build_product_interface: iface.glob2hat must be integer-valued.');
    end
    if any(g2h(:) < 0)
      error('build_product_interface: iface.glob2hat must be nonnegative (0..nHat).');
    end
    if any(g2h(:) > nHat)
      error('build_product_interface: iface.glob2hat contains values > iface.nHat.');
    end
  end
  g2h = g2h(:); % treat as column
  mapLen = numel(g2h);

  % CHANGED: if nHat==0, mapping must not claim any assembled DOF
  if nHat == 0 && ~isempty(g2h) && any(g2h(:) ~= 0)
    error('build_product_interface: iface.nHat==0 but iface.glob2hat contains nonzero entries.');
  end

  nSub = numel(sub);

  % Count total number of product DOFs (i.e. product space dimension)
  nProd = 0;
  for i = 1:nSub
    if ~isfield(sub(i),'glob_G')
      error('build_product_interface: sub(%d) missing glob_G (run identify_interface_dofs).', i);
    end

    % CHANGED: validate sub(i).glob_G upfront (type + shape), reject char/logical
    gG = sub(i).glob_G;
    if ischar(gG) || islogical(gG) || (exist('isstring','file')==2 && isstring(gG))
      error('build_product_interface: sub(%d).glob_G must be numeric (char/logical/string not allowed).', i);
    end
    if ~(isempty(gG) || (isnumeric(gG) && isreal(gG) && isvector(gG) && all(isfinite(gG(:)))))
      error('build_product_interface: sub(%d).glob_G must be a real, finite numeric vector (or empty).', i);
    end
    if ~isempty(gG)
      if any(gG(:) ~= fix(gG(:)))
        error('build_product_interface: sub(%d).glob_G must be integer-valued.', i);
      end
      if any(gG(:) < 1)
        error('build_product_interface: sub(%d).glob_G must be in 1..numel(iface.glob2hat).', i);
      end
      if mapLen < 1
        error('build_product_interface: iface.glob2hat is empty but sub(%d).glob_G is nonempty.', i);
      end
      if any(gG(:) > mapLen)
        error('build_product_interface: sub(%d).glob_G contains out-of-range indices for iface.glob2hat.', i);
      end
      if numel(unique(gG(:))) ~= numel(gG(:))
        error('build_product_interface: sub(%d).glob_G contains duplicates (ambiguous indexing).', i);
      end
    end

    nProd = nProd + numel(gG);
  end

  prod = struct();
  prod.nHat = nHat;
  prod.nProd = nProd;
  prod.prod2sub = zeros(nProd, 1);
  prod.prod2hat = zeros(nProd, 1);
  prod.hat2prod = cell(nHat, 1);

  % Assign product indices by concatenating subdomains 1..nSub.
  cursor = 1;
  for i = 1:nSub
    gG = sub(i).glob_G(:); % global free DOFs on interface for subdomain i (column)

    if isempty(gG)
      sub(i).gamma_glob = gG;
      sub(i).gamma_hat  = zeros(0,1);
      sub(i).prod_idx   = zeros(0,1);
      continue;
    end

    % CHANGED: robust mapping with additional sanity checks
    hat = g2h(gG); % map global free DOFs -> assembled interface indices

    if any(hat == 0)
      error('build_product_interface: sub(%d) contains a glob_G DOF that is not interface (glob2hat==0).', i);
    end
    if any(hat < 1) || any(hat > nHat)
      error('build_product_interface: sub(%d) maps to invalid hat indices (outside 1..nHat).', i);
    end
    if numel(unique(hat)) ~= numel(hat)
      % CHANGED: disallow duplicate hat mapping within one subdomain (typically caused by bad glob_G)
      error('build_product_interface: sub(%d) maps multiple glob_G entries to the same hat index (duplicate hat).', i);
    end

    idx = (cursor:(cursor+numel(gG)-1)).';
    cursor = cursor + numel(gG);

    sub(i).gamma_glob = gG;
    sub(i).gamma_hat  = hat(:);
    sub(i).prod_idx   = idx;

    prod.prod2sub(idx) = i;
    prod.prod2hat(idx) = hat(:);

    % Fill hat2prod lists (inverse mapping: hat DOF -> product copies)
    for k = 1:numel(hat)
      hk = hat(k);
      prod.hat2prod{hk}(end+1,1) = idx(k); %#ok<AGROW>
    end
  end

  % Basic sanity.
  if cursor ~= nProd + 1
    error('build_product_interface: internal error in product indexing.');
  end
end