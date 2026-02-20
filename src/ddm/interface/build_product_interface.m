function [sub, prod] = build_product_interface(sub, iface)
%BUILD_PRODUCT_INTERFACE  Construct product (duplicated) interface indexing and bookkeeping.
%
% Link to thesis:
%   Chapter 3.3.2 (product interface space), equations (3.45)–(3.47).
%
% This routine constructs the product interface space W = Π_i W_i by assigning
% a unique "product" index to every local interface DOF on every subdomain.
% Physical interface DOFs (assembled interface space \hat{W}) are represented
% by iface.hat2glob / iface.glob2hat.
%
% Inputs:
%   sub   : subdomain array with fields dofs_G (local indices) and glob_G (global free DOFs).
%   iface : interface bookkeeping from identify_interface_dofs (needs glob2hat and nHat).
%
% Outputs:
%   sub  : updated with fields:
%     .gamma_glob : local interface DOFs as global free DOFs (same as .glob_G)
%     .gamma_hat  : corresponding assembled interface indices (1..nHat)
%     .prod_idx   : product indices for this subdomain's interface vector (contiguous segment)
%
%   prod : product-space bookkeeping:
%     .nHat      : number of assembled interface DOFs
%     .nProd     : dimension of product space (sum_i |F_G^(i)|)
%     .prod2sub  : (nProd x 1) subdomain id for each product DOF
%     .prod2hat  : (nProd x 1) assembled interface index for each product DOF
%     .hat2prod  : cell (nHat x 1): hat2prod{k} lists product copies of hat DOF k

  if nargin ~= 2
    error('build_product_interface: expected inputs (sub, iface).');
  end
  if ~isstruct(iface) || ~isfield(iface,'glob2hat') || ~isfield(iface,'nHat')
    error('build_product_interface: iface must contain glob2hat and nHat (run identify_interface_dofs).');
  end

  nSub = numel(sub);
  nHat = iface.nHat;

  % Count total number of product DOFs (i.e. compute product space dimension)
  nProd = 0;
  for i = 1:nSub
    if ~isfield(sub(i),'glob_G')
      error('build_product_interface: sub(%d) missing glob_G (run identify_interface_dofs).', i);
    end
    nProd = nProd + numel(sub(i).glob_G);
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
    gG = sub(i).glob_G(:);                 % global free DOFs on interface for subdomain i
    if isempty(gG)
      sub(i).gamma_glob = gG;
      sub(i).gamma_hat  = zeros(0,1);
      sub(i).prod_idx   = zeros(0,1);
      continue;
    end

    hat = iface.glob2hat(gG); %map global free DOFs on interface for subdomain i to hat indices
    if any(hat == 0)
      error('build_product_interface: sub(%d) contains a glob_G DOF that is not recognized as interface in iface.glob2hat.', i);
    end

    idx = (cursor:(cursor+numel(gG)-1)).';
    cursor = cursor + numel(gG);

    sub(i).gamma_glob = gG;
    sub(i).gamma_hat  = hat;
    sub(i).prod_idx   = idx;

    prod.prod2sub(idx) = i;
    prod.prod2hat(idx) = hat;

    % Fill hat2prod lists.
    % the inverse mapping: hat DOF → all product copies that represent it
    % for each assembled interface DOF hk we get product indices of every duplicated copy of that DOF across all adjacent subdomains
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
