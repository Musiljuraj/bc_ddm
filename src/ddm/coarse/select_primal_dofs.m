function primal = select_primal_dofs(p, ddm, iface)
%SELECT_PRIMAL_DOFS  Select primal (coarse) interface DOFs for a structured subdomain decomposition.
%
% Link to thesis:
%   Chapter 3.4.2 (choice of primal variables: corner nodes of the decomposition),
%   equations (3.58)–(3.60).
%
% In this thesis' model setting (unit square, rectangular subdomains), primal DOFs
% are chosen as "corner nodes of the decomposition", i.e. interface nodes whose
% coordinates lie on both a vertical and a horizontal subdomain partition line:
%   x = k / nSubX,  y = l / nSubY
% for integers k,l.
%
% Inputs:
%   p    : Np x 2 node coordinates.
%   ddm  : struct from build_subdomains_structured (needs nSubX,nSubY,dof2node).
%   iface: struct from identify_interface_dofs (needs glob and glob2hat).
%
% Output:
%   primal : struct with fields:
%       .glob_c     : global free DOFs selected as primal
%       .glob_delta : remaining interface global free DOFs
%       .hat_c      : assembled interface indices of primal DOFs
%       .hat_delta  : assembled interface indices of remaining interface DOFs
%       .nC         : number of primal DOFs

  if nargin ~= 3
    error('select_primal_dofs: expected inputs (p, ddm, iface).');
  end

  % Fixed coordinate tolerance (no opts input).
  tol = 1e-12;

  if ~isfield(ddm,'nSubX') || ~isfield(ddm,'nSubY') || ~isfield(ddm,'dof2node')
    error('select_primal_dofs: ddm must contain nSubX, nSubY, dof2node.');
  end
  if ~isfield(iface,'glob') || ~isfield(iface,'glob2hat')
    error('select_primal_dofs: iface must contain glob and glob2hat.');
  end

  nSubX = ddm.nSubX;
  nSubY = ddm.nSubY;

  %lines that defines subdomain boundaries (their cross-points are corners)
  xLines = (0:nSubX) / nSubX;
  yLines = (0:nSubY) / nSubY;

  iface_glob = iface.glob(:);
  nHat = numel(iface_glob);

  is_primal = false(nHat, 1); 

  %for all dofs in w-hat find out whether they lay on a cross-point of subdomain boundary lines
  for k = 1:nHat
    g = iface_glob(k);                 % global free DOF id
    node = ddm.dof2node(g);            % global node index
    x = p(node,1); 
    y = p(node,2);

    on_x = any(abs(x - xLines) < tol);
    on_y = any(abs(y - yLines) < tol);

    if on_x && on_y
      is_primal(k) = true;
    end
  end

  
  primal = struct();
  primal.glob_c = sort(iface_glob(is_primal));
  primal.glob_delta = sort(iface_glob(~is_primal));

  primal.hat_c = iface.glob2hat(primal.glob_c);
  primal.hat_delta = iface.glob2hat(primal.glob_delta);

  primal.nC = numel(primal.glob_c);

  % Store for convenience.
  primal.tol = tol;
end