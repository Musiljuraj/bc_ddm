function sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f_handle)
%ASSEMBLE_SUBDOMAIN_MATRICES_P1  Assemble local subdomain stiffness matrices and load vectors (P1).
%
% Link to thesis:
%   Chapter 3.2.1 (local subdomain systems), equations (3.11)–(3.16).
%
% For each subdomain Ω_i (represented by the element index set sub(i).elems),
% we assemble a local stiffness matrix K^(i) and a local load vector f^(i) by
% summing element contributions belonging to Ω_i only.
%
% IMPORTANT (theory-consistent):
%   K^(i) must be assembled from the elements in Ω_i, not obtained by slicing
%   a globally assembled K. This matches the thesis definition and avoids
%   double-counting interface contributions.
%
% Dirichlet handling:
%   Local matrices/vectors are assembled only for *free* DOFs (Dirichlet nodes
%   are skipped using ddm.node2dof==0).
%
% Inputs:
%   p        : Np x 2 node coordinates.
%   t        : Nt x 3 connectivity (global node indices).
%   sub      : subdomain array from build_subdomains_structured (needs .elems, .glob2loc, .nloc).
%   ddm      : global bookkeeping (needs .node2dof).
%   f_handle : function handle f(x,y). If omitted, uses f=0.
%
% Output:
%   sub : updated with fields:
%       .K : nloc x nloc sparse stiffness matrix
%       .f : nloc x 1 sparse load vector

  if nargin < 4 || nargin > 5
    error('assemble_subdomain_matrices_P1: expected inputs (p,t,sub,ddm[,f_handle]).');
  end
  if nargin == 4
    f_handle = @(x,y) 0;
  end
  if ~isa(f_handle, 'function_handle')
    error('assemble_subdomain_matrices_P1: f_handle must be a function handle.');
  end
  if ~isstruct(ddm) || ~isfield(ddm,'node2dof')
    error('assemble_subdomain_matrices_P1: ddm must contain field node2dof.');
  end

  nSub = numel(sub);

  for i = 1:nSub
    nloc = sub(i).nloc;

    % Trivial case: no free DOFs in this subdomain.
    if nloc == 0
      sub(i).K = sparse(0,0);
      sub(i).f = sparse(0,1);
      continue;
    end

    elems = sub(i).elems(:);
    nE = numel(elems);

    % Preallocate triplets  for local K_e^(i)(upper bound; some entries skipped due to Dirichlet vertices).
    I = zeros(9*nE, 1);
    J = zeros(9*nE, 1);
    V = zeros(9*nE, 1);
    idxK = 1;
    %preallocate for right hand side (local element load vector)
    IF = zeros(3*nE, 1);
    VF = zeros(3*nE, 1);
    idxF = 1;

    for ee = 1:nE
      e = elems(ee);
      nodes = t(e,:);          % 1x3 global node indices
      xy = p(nodes,:);         % 3x2 coordinates

      Ke = triP1_stiffness(xy);
      fe = triP1_load(xy, f_handle);

      dofs = ddm.node2dof(nodes);  % 1x3 global free DOF ids (0 for Dirichlet)

      % Stiffness load contribution: local indices (0 if Dirichlet or not in local set - should not happen for free nodes in elements)
      loc = zeros(1,3);
      for a = 1:3
        if dofs(a) > 0
          la = sub(i).glob2loc(dofs(a)); %local index of particular node a in subdomain right hand side vector
          if la <= 0
            error('assemble_subdomain_matrices_P1: sub(%d) is missing a free DOF that appears in its element set (global dof=%d).', i, dofs(a));
          end
          loc(a) = la;

          IF(idxF) = la;
          VF(idxF) = fe(a);
          idxF = idxF + 1;
        end
      end

      % Stiffness K contributions (skip Dirichlet vertices).
      for a = 1:3
        if dofs(a) == 0, continue; end
        la = loc(a); %local row index of particular node a in subdomain Ke
        for b = 1:3
          if dofs(b) == 0, continue; end
          lb = loc(b); %local column index of particular node a in subdomain Ke

          I(idxK) = la;
          J(idxK) = lb;
          V(idxK) = Ke(a,b);
          idxK = idxK + 1;
        end
      end
    end

    % Trim triplets to used length.
    I = I(1:idxK-1); J = J(1:idxK-1); V = V(1:idxK-1);
    IF = IF(1:idxF-1); VF = VF(1:idxF-1);

    sub(i).K = sparse(I, J, V, nloc, nloc);
    sub(i).f = sparse(IF, ones(size(IF)), VF, nloc, 1);
  end
end
