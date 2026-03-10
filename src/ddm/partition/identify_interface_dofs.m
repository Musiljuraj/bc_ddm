function [sub, iface] = identify_interface_dofs(sub, ddm)
%IDENTIFY_INTERFACE_DOFS Identify interior and interface DOFs of each subdomain.
% Thesis link: Chapter 4.1.3 (interior/interface splitting).
% This routine classifies local unknowns for later Schur complement reduction.
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
%       .nI, .nG : counts of interior/interface DOFs
%
%   iface : global interface bookkeeping:
%       .counts    : nFree x 1 multiplicity of each global free DOF across subdomains
%       .glob      : list of global free DOFs that are shared (counts>1)
%       .nHat      : number of unique (assembled) interface DOFs
%       .glob2hat  : mapping global -> hat indexing in assembled interface space \hat{W}
%       .hat2glob  : mapping hat -> global (equals iface.glob)
%       .dof2sub   : cell array; for each interface DOF (hat index), list of subdomains containing it

  if nargin ~= 2
    error('identify_interface_dofs: expected inputs (sub, ddm).');
  end

  % FIXED: validate 'sub' type early (prevents non-struct/odd inputs)
  if ~isstruct(sub) % FIXED
    error('identify_interface_dofs: sub must be a struct array.'); % FIXED
  end

  % FIXED: validate ddm.nFree existence and contract (scalar positive integer)
  if ~isstruct(ddm) || ~isfield(ddm,'nFree')
    error('identify_interface_dofs: ddm must contain field nFree.');
  end
  nFree = ddm.nFree; % CHANGED
  if ~(isnumeric(nFree) && isreal(nFree) && isfinite(nFree) && isscalar(nFree) && nFree >= 1 && nFree == fix(nFree)) % FIXED
    error('identify_interface_dofs: ddm.nFree must be a positive scalar integer.'); % FIXED
  end
  nFree = double(nFree); % ADDED (robust for integer class)

  nSub = numel(sub);

  % ADDED: pre-validate all sub(i).loc2glob to reject char/logical/NaN/Inf/complex and duplicates
  for i = 1:nSub
    if ~isfield(sub(i),'loc2glob') % ADDED
      error('identify_interface_dofs: sub(%d) must contain field loc2glob.', i); % ADDED
    end

    gi_raw = sub(i).loc2glob; % ADDED

    % FIXED: reject char inputs (ASCII indices) and logical inputs (0/1 indices)
    if ~(isnumeric(gi_raw) && ~islogical(gi_raw)) % FIXED
      error('identify_interface_dofs: sub(%d).loc2glob must be numeric (not char/logical).', i); % FIXED
    end
    if ~isreal(gi_raw) % ADDED
      error('identify_interface_dofs: sub(%d).loc2glob must be real-valued.', i); % ADDED
    end
    if ~isvector(gi_raw) % ADDED
      error('identify_interface_dofs: sub(%d).loc2glob must be a vector.', i); % ADDED
    end

    gi = gi_raw(:); % ADDED
    if issparse(gi) % ADDED
      gi = full(gi); % ADDED
    end

    % FIXED: reject NaN/Inf explicitly (otherwise comparisons may miss NaN)
    if any(~isfinite(gi)) % FIXED
      error('identify_interface_dofs: sub(%d).loc2glob must be finite.', i); % FIXED
    end

    % FIXED: enforce integer-valued indices
    if any(gi ~= fix(gi)) % FIXED
      error('identify_interface_dofs: sub(%d).loc2glob must be integer-valued.', i); % FIXED
    end

    % Existing range check, but now guaranteed meaningful
    if any(gi < 1) || any(gi > nFree)
      error('identify_interface_dofs: sub(%d).loc2glob contains invalid global DOF indices.', i);
    end

    % FIXED: guard against duplicates within a subdomain (would inflate multiplicities and dof2sub)
    if numel(unique(gi)) ~= numel(gi) % FIXED
      error('identify_interface_dofs: sub(%d).loc2glob contains duplicate global DOF indices.', i); % FIXED
    end

    % ADDED: normalize stored loc2glob to a column vector for downstream consistency
    sub(i).loc2glob = gi; % ADDED
  end

  % -----------------------------
  % 1) Multiplicity of each global free DOF across subdomains
  % -----------------------------
  counts = zeros(nFree, 1);
  for i = 1:nSub
    gi = sub(i).loc2glob; % CHANGED (already validated/normalized)
    counts(gi) = counts(gi) + 1;
  end

  % Build iface structure
  iface = struct();
  iface.counts = counts;
  iface.glob = find(counts > 1); % global interface dofs (shared)
  iface.nHat = numel(iface.glob);

  % Mapping hat<->glob
  glob2hat = zeros(nFree, 1);
  glob2hat(iface.glob) = (1:iface.nHat).';
  iface.glob2hat = glob2hat;
  iface.hat2glob = iface.glob;

  % Build dof2sub lists.
  iface.dof2sub = cell(iface.nHat, 1);
  for i = 1:nSub
    gi = sub(i).loc2glob; % ADDED
    maskG = counts(gi) > 1; % ADDED
    gG = gi(maskG);         % CHANGED (clearer than indexing sub(i).loc2glob twice)
    gG = unique(gG);        % FIXED (defensive: prevents duplicates in dof2sub)
    hat_ids = glob2hat(gG);
    for k = 1:numel(hat_ids)
      iface.dof2sub{hat_ids(k)}(end+1,1) = i; %#ok<AGROW>
    end
  end

  % FIXED: ensure dof2sub entries are unique (defensive) and consistent length if multiplicities are correct
  for k = 1:iface.nHat % FIXED
    if ~isempty(iface.dof2sub{k}) % FIXED
      iface.dof2sub{k} = unique(iface.dof2sub{k}(:), 'stable'); % FIXED
    end
  end

  % -----------------------------
  % 2) Per-subdomain split into interior vs interface
  % -----------------------------
  for i = 1:nSub
    gi = sub(i).loc2glob;      % CHANGED (validated/normalized)
    maskG = counts(gi) > 1;

    sub(i).dofs_G = find(maskG);
    sub(i).dofs_I = find(~maskG);

    sub(i).glob_G = gi(maskG);
    sub(i).glob_I = gi(~maskG);

    sub(i).nG = numel(sub(i).dofs_G);
    sub(i).nI = numel(sub(i).dofs_I);
  end
end