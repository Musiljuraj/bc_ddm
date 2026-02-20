function [sub, iface] = identify_interface_dofs(sub, ddm)
%IDENTIFY_INTERFACE_DOFS  Identify interior/interface DOFs per subdomain and global interface bookkeeping.
%
% Link to thesis:
%   Chapter 3.1.3 (interior vs interface DOFs), equation (3.10).
%
% Discrete interpretation used here (robust for FE meshes):
%   - A *free* global DOF is considered an interface DOF if it appears in the
%     local DOF sets of at least two subdomains (i.e., it is shared).
%   - All other free DOFs are treated as "interior" with respect to the
%     interface reduction (they do not couple multiple subdomains).
%
% Inputs:
%   sub : subdomain struct array from build_subdomains_structured (must contain .loc2glob).
%   ddm : global bookkeeping struct (must contain .nFree).
%
% Outputs:
%   sub : updated with per-subdomain splits:
%       .dofs_I  : local indices of interior DOFs in the full local vector
%       .dofs_G  : local indices of interface DOFs in the full local vector
%       .glob_I  : corresponding global free DOFs
%       .glob_G  : corresponding global free DOFs
%       .nI, .nG : counts of interface and inner nodes
%
%   iface : global interface bookkeeping:
%       .counts    : nFree x 1 multiplicity of each global free DOF across subdomains
%       .glob      : list of global free DOFs that are shared (counts>1)
%       .nHat      : number of unique (assembled) interface DOFs
%       .glob2hat  : mapping of nodes from global indexing to assembled interface space overhat{W} indexing
%       .hat2glob  : mapping of nodes from assembled interface space overhat{W} to global indexing 
%       .dof2sub   : cell array, lists of indices of a subdomains, to which each interface node belong

  if nargin ~= 2
    error('identify_interface_dofs: expected inputs (sub, ddm).');
  end
  if ~isstruct(ddm) || ~isfield(ddm,'nFree')
    error('identify_interface_dofs: ddm must contain field nFree.');
  end

  nSub = numel(sub);
  nFree = ddm.nFree;

  % -----------------------------
  % 1) Multiplicity of each global free DOF across subdomains
  % check each node occurence in every subdomain, for each occurence increment <counts> by one at the position of this node
  % -----------------------------
  counts = zeros(nFree, 1);
  for i = 1:nSub
    gi = sub(i).loc2glob(:);
    if any(gi < 1) || any(gi > nFree)
      error('identify_interface_dofs: sub(%d).loc2glob contains invalid global DOF indices.', i);
    end
    counts(gi) = counts(gi) + 1;
  end

  %Build iface structure 
  iface = struct();
  iface.counts = counts; 
  iface.glob = find(counts > 1); %dof's (in global indx) of the assembled interface space overhat{W} 
  iface.nHat = numel(iface.glob);  

  %mapping of nodes from assembled interface space overhat{W} to global indexing and (hat2glob) and inverse mapping (glob2hat)
  glob2hat = zeros(nFree,1);
  glob2hat(iface.glob) = (1:iface.nHat).';
  iface.glob2hat = glob2hat;
  iface.hat2glob = iface.glob;

  % Build dof2sub lists.
  % Create list of variable-lenght lists each containing indices of a subdomains, to which each interface node belong
  iface.dof2sub = cell(iface.nHat, 1);
  for i = 1:nSub
    gG = sub(i).loc2glob(counts(sub(i).loc2glob) > 1); %global indexes of interface DOFs in subdomain i
    hat_ids = glob2hat(gG); %mapp these indices into assembled interface space indices
    for k = 1:numel(hat_ids)
      iface.dof2sub{hat_ids(k)}(end+1,1) = i; %adding number of domian to which each k node belong #ok<AGROW>
    end
  end

  % -----------------------------
  % 2) Per-subdomain split into interior vs interface
  % -----------------------------
  for i = 1:nSub
    gi = sub(i).loc2glob(:); %global indices of nodes in subdomain
    maskG = counts(gi) > 1; % mask of those which belongs to multiple subdomains

    sub(i).dofs_G = find(maskG); 
    sub(i).dofs_I = find(~maskG); 

    sub(i).glob_G = gi(maskG); %interface nodes of subdomain i, in global indexing
    sub(i).glob_I = gi(~maskG); %inner nodes of subdomain i, in global indexing

    sub(i).nG = numel(sub(i).dofs_G);
    sub(i).nI = numel(sub(i).dofs_I);
  end
end
