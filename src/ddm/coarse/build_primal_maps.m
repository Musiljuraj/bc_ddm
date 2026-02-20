function [sub, primal] = build_primal_maps(sub, ddm, iface, prod, primal)
%BUILD_PRIMAL_MAPS  Build per-subdomain primal/delta index sets and global primal indexing.
%
% Link to thesis:
%   Chapter 3.4 (primal variables and interface split), equation (3.58),
%   and the decomposition W = W^c ⊕ W^Δ (3.61) with reduced constraints (3.63).
%
% Inputs:
%   sub   : subdomain array with fields:
%           .glob_G    : global free DOFs on interface (local order)
%           .prod_idx  : product indices (from build_product_interface)
%   ddm   : global bookkeeping (needs .nFree)
%   iface : interface bookkeeping (needs .glob2hat and .nHat)
%   prod  : product bookkeeping (needs .nProd, .prod2hat)
%   primal: primal selection (from select_primal_dofs), needs .glob_c and .hat_c
%
% Outputs:
%   sub : updated with per-subdomain sets:
%       .idx_c         : indices of primal DOFs in the local interface vector (1..nG)
%       .idx_d         : indices of delta DOFs in the local interface vector (1..nG)
%       .glob_c        : global free DOFs of primal interface nodes on this subdomain
%       .glob_d        : global free DOFs of delta interface nodes on this subdomain
%       .c_ids         : global primal ids (1..nC) corresponding to glob_c
%       .prod_idx_c    : product indices corresponding to primal interface DOFs
%       .prod_idx_d    : product indices corresponding to delta interface DOFs
%
%   primal : enriched with global maps:
%       .glob2c        : (nFree x 1) map global free DOF -> primal id (0 if not primal)
%       .hat2c         : (nHat x 1)  map assembled interface index -> primal id (0 if not primal)
%       .prod_is_c     : (nProd x 1) logical, true if product DOF corresponds to a primal physical DOF
%       .prod_is_d     : (nProd x 1) logical complement
%       .delta_idx     : list of product indices belonging to delta part
%       .prod2delta    : (nProd x 1) map product index -> delta index (0 if primal)
%       .nDeltaProd    : number of delta product DOFs

  if nargin ~= 5
    error('build_primal_maps: expected inputs (sub, ddm, iface, prod, primal).');
  end
  if ~isfield(ddm,'nFree')
    error('build_primal_maps: ddm must contain nFree.');
  end
  if ~isfield(iface,'nHat') || ~isfield(iface,'glob2hat')
    error('build_primal_maps: iface must contain nHat and glob2hat.');
  end
  if ~isfield(prod,'nProd') || ~isfield(prod,'prod2hat')
    error('build_primal_maps: prod must contain nProd and prod2hat.');
  end
  if ~isfield(primal,'glob_c')
    error('build_primal_maps: primal must contain glob_c (run select_primal_dofs).');
  end

  %sizes of all three index spaces (free DOFs, W-hat and W-prod)
  nFree = ddm.nFree;
  nHat  = iface.nHat;
  nProd = prod.nProd;

  % Global primal indexing.
  glob_c = primal.glob_c(:);
  nC = numel(glob_c);

  glob2c = zeros(nFree, 1);
  if nC > 0
    glob2c(glob_c) = (1:nC).';
  end
  primal.glob2c = glob2c; %primal indexing in global range (i.e. 0 if not primal, 1 .. nC if primal)
  primal.nC = nC;

  % Map assembled interface DOF -> primal id.
  hat2c = zeros(nHat, 1);
  if isfield(primal,'hat_c') && ~isempty(primal.hat_c)
    hat2c(primal.hat_c(:)) = (1:nC).';
  else
    % In case hat_c not present, compute from glob_c.
    hat2c(iface.glob2hat(glob_c)) = (1:nC).';
  end
  primal.hat2c = hat2c; %primal indexing in W-hat range (i.e. 0 if not primal, 1 .. nC if primal)

  % Per-subdomain split.
  nSub = numel(sub);
  prod_is_c = false(nProd, 1);

  for i = 1:nSub
    if ~isfield(sub(i),'gamma_glob')
      % gamma_glob is identical to glob_G; allow either field name.
      if isfield(sub(i),'glob_G')
        sub(i).gamma_glob = sub(i).glob_G(:);
      else
        error('build_primal_maps: sub(%d) must contain glob_G or gamma_glob.', i);
      end
    end
    if ~isfield(sub(i),'prod_idx')
      error('build_primal_maps: sub(%d) missing prod_idx (run build_product_interface).', i);
    end

    gG = sub(i).gamma_glob(:);
    nG = numel(gG);

    if nG == 0
      sub(i).idx_c = zeros(0,1);
      sub(i).idx_d = zeros(0,1);
      sub(i).glob_c = zeros(0,1);
      sub(i).glob_d = zeros(0,1);
      sub(i).c_ids  = zeros(0,1);
      sub(i).prod_idx_c = zeros(0,1);
      sub(i).prod_idx_d = zeros(0,1);
      continue;
    end

    c_ids_local = glob2c(gG);          % 0 if not primal
    is_c = (c_ids_local > 0);           %mask for local's 

    sub(i).idx_c = find(is_c); % Store local positions of primal entries in the local interface vector ordering.
    sub(i).idx_d = find(~is_c); % Store local positions of delta entries in the local interface vector ordering.

    sub(i).glob_c = gG(is_c); % Store global positions of primal entries in the local interface vector ordering.
    sub(i).glob_d = gG(~is_c); % Store global positions of delta entries in the local interface vector ordering.

    sub(i).c_ids  = c_ids_local(is_c);  % Store the coarse ids corresponding to the subdomain’s primal interface DOFs,

    % Split product indices accordingly.
    pidx = sub(i).prod_idx(:);
    if numel(pidx) ~= nG
      error('build_primal_maps: sub(%d) prod_idx length mismatch with gamma_glob.', i);
    end
    sub(i).prod_idx_c = pidx(is_c);
    sub(i).prod_idx_d = pidx(~is_c);

    prod_is_c(sub(i).prod_idx_c) = true;
  end

  primal.prod_is_c = prod_is_c;
  primal.prod_is_d = ~prod_is_c;

  delta_idx = find(~prod_is_c);
  primal.delta_idx = delta_idx;
  primal.nDeltaProd = numel(delta_idx);

  prod2delta = zeros(nProd, 1);
  prod2delta(delta_idx) = (1:primal.nDeltaProd).';
  primal.prod2delta = prod2delta;
end
