function sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f_handle)
%ASSEMBLE_SUBDOMAIN_MATRICES_P1 Assemble local FEM matrices on each subdomain.
% Thesis link: Chapter 4.1.2 and 4.2.1 (local spaces and local systems).
% The routine forms subdomain-wise stiffness and load data induced by the
% global P1 discretization.

  % -------------------- signature / basic checks --------------------
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
  if ~isstruct(sub)
    error('assemble_subdomain_matrices_P1: sub must be a struct array.');
  end

  % -------------------- input hardening (mesh + ddm) --------------------
  % Validate p
  if ~isnumeric(p) || ndims(p) ~= 2 || size(p,2) ~= 2
    error('assemble_subdomain_matrices_P1: p must be an Np-by-2 numeric matrix.');
  end
  if ~isreal(p) || any(~isfinite(p(:)))
    error('assemble_subdomain_matrices_P1: p must be real and finite.');
  end
  Np = size(p,1);

  % Validate t
  if ~isnumeric(t) || ndims(t) ~= 2 || size(t,2) ~= 3
    error('assemble_subdomain_matrices_P1: t must be an Nt-by-3 numeric matrix.');
  end
  if ~isreal(t) || any(~isfinite(t(:)))
    error('assemble_subdomain_matrices_P1: t must be real and finite.');
  end
  if any(t(:) ~= fix(t(:)))
    error('assemble_subdomain_matrices_P1: t must contain integer-valued node indices.');
  end
  if any(t(:) < 1) || any(t(:) > Np)
    error('assemble_subdomain_matrices_P1: t contains out-of-range node indices (must be in 1..size(p,1)).');
  end
  Nt = size(t,1);

  % Validate ddm.node2dof
  node2dof = ddm.node2dof(:);
  if ~isnumeric(node2dof) || ~isreal(node2dof) || any(~isfinite(node2dof))
    error('assemble_subdomain_matrices_P1: ddm.node2dof must be numeric, real, finite.');
  end
  if numel(node2dof) ~= Np
    error('assemble_subdomain_matrices_P1: ddm.node2dof must have length size(p,1).');
  end
  if any(node2dof ~= fix(node2dof)) || any(node2dof < 0)
    error('assemble_subdomain_matrices_P1: ddm.node2dof must be integer-valued and nonnegative (0 for Dirichlet).');
  end

  if isfield(ddm,'nFree') && isnumeric(ddm.nFree) && isscalar(ddm.nFree) && isfinite(ddm.nFree)
    nFree = ddm.nFree;
    if ~isreal(nFree) || nFree ~= fix(nFree) || nFree < 0
      error('assemble_subdomain_matrices_P1: ddm.nFree must be a nonnegative integer scalar.');
    end
  else
    nFree = max(node2dof);
  end
  if any(node2dof > nFree)
    error('assemble_subdomain_matrices_P1: ddm.node2dof contains values larger than nFree.');
  end

  % Light f_handle sanity: scalar inputs should return a real finite scalar numeric.
  try
    tmp = f_handle(p(1,1), p(1,2));
  catch
    error('assemble_subdomain_matrices_P1: f_handle must accept scalar (x,y) inputs.');
  end
  if ~isnumeric(tmp) || ~isscalar(tmp) || ~isreal(tmp) || ~isfinite(tmp)
    error('assemble_subdomain_matrices_P1: f_handle(x,y) must return a real finite scalar numeric value for scalar inputs.');
  end

  % -------------------- subdomain contract (relaxed glob2loc invariant) --------------------
  nSub = numel(sub);
  for i = 1:nSub
    if ~isfield(sub(i),'nloc') || ~isfield(sub(i),'elems') || ~isfield(sub(i),'glob2loc')
      error('assemble_subdomain_matrices_P1: sub(%d) must contain fields nloc, elems, glob2loc.', i);
    end

    nloc = sub(i).nloc;
    if ~isnumeric(nloc) || ~isscalar(nloc) || ~isfinite(nloc) || ~isreal(nloc) || nloc ~= fix(nloc) || nloc < 0
      error('assemble_subdomain_matrices_P1: sub(%d).nloc must be a nonnegative integer scalar.', i);
    end

    elems = sub(i).elems(:);
    if ~isnumeric(elems) || ~isreal(elems) || any(~isfinite(elems)) || any(elems ~= fix(elems))
      error('assemble_subdomain_matrices_P1: sub(%d).elems must be integer-valued, real, finite.', i);
    end
    if any(elems < 1) || any(elems > Nt)
      error('assemble_subdomain_matrices_P1: sub(%d).elems contains out-of-range element indices (must be in 1..size(t,1)).', i);
    end
    if numel(elems) ~= numel(unique(elems))
      error('assemble_subdomain_matrices_P1: sub(%d).elems contains duplicates (would double-count element contributions).', i);
    end
    if nloc > 0 && isempty(elems)
      error('assemble_subdomain_matrices_P1: sub(%d) has nloc>0 but empty elems (inconsistent subdomain definition).', i);
    end

    g2l = sub(i).glob2loc(:);
    if ~isnumeric(g2l) || ~isvector(g2l) || ~isreal(g2l) || any(~isfinite(nonzeros(g2l)))
      error('assemble_subdomain_matrices_P1: sub(%d).glob2loc must be a real finite numeric vector (zeros allowed).', i);
    end

    % Check mapping only for DOFs that actually appear in this subdomain’s element set.
    if ~isempty(elems)
      nodes_used = unique(t(elems, :)(:));
      dofs_used = unique(node2dof(nodes_used));
      dofs_used = dofs_used(dofs_used > 0); % free DOFs only

      if ~isempty(dofs_used)
        if numel(g2l) < max(dofs_used)
          error('assemble_subdomain_matrices_P1: sub(%d).glob2loc is too short to index used DOFs from elems.', i);
        end

        la = g2l(dofs_used);

        if any(la <= 0)
          error('assemble_subdomain_matrices_P1: sub(%d).glob2loc is missing a free DOF used by elems.', i);
        end
        if any(la ~= fix(la))
          error('assemble_subdomain_matrices_P1: sub(%d).glob2loc must map used DOFs to integer local indices.', i);
        end
        if any(la < 1) || any(la > nloc)
          error('assemble_subdomain_matrices_P1: sub(%d).glob2loc maps a used DOF to an index outside 1..nloc.', i);
        end
        if numel(unique(la)) ~= numel(la)
          error('assemble_subdomain_matrices_P1: sub(%d).glob2loc maps distinct used DOFs to the same local index.', i);
        end
      else
        % All nodes in elems are Dirichlet (rare). Then assembly should produce 0x0 or consistent empties.
        if nloc ~= 0
          error('assemble_subdomain_matrices_P1: sub(%d) elems contain no free DOFs but nloc>0 (inconsistent).', i);
        end
      end
    end

    % Optional: if loc2glob exists, check basic sanity (non-fatal for assembly, but useful)
    if isfield(sub(i),'loc2glob') && ~isempty(sub(i).loc2glob)
      l2g = sub(i).loc2glob(:);
      if ~isnumeric(l2g) || ~isreal(l2g) || any(~isfinite(l2g)) || any(l2g ~= fix(l2g))
        error('assemble_subdomain_matrices_P1: sub(%d).loc2glob must be integer-valued, real, finite.', i);
      end
      if numel(l2g) ~= nloc
        error('assemble_subdomain_matrices_P1: sub(%d).loc2glob length must equal nloc.', i);
      end
      if any(l2g < 1) || any(l2g > nFree)
        error('assemble_subdomain_matrices_P1: sub(%d).loc2glob contains out-of-range global DOFs (must be in 1..nFree).', i);
      end
      if numel(unique(l2g)) ~= numel(l2g)
        error('assemble_subdomain_matrices_P1: sub(%d).loc2glob contains duplicates (invalid mapping).', i);
      end
    end
  end

  % -------------------- assembly --------------------
  for i = 1:nSub
    nloc = sub(i).nloc;

    if nloc == 0
      sub(i).K = sparse(0,0);
      sub(i).f = sparse(0,1);
      continue;
    end

    elems = sub(i).elems(:);
    nE = numel(elems);

    I = zeros(9*nE, 1);
    J = zeros(9*nE, 1);
    V = zeros(9*nE, 1);
    idxK = 1;

    IF = zeros(3*nE, 1);
    VF = zeros(3*nE, 1);
    idxF = 1;

    for ee = 1:nE
      e = elems(ee);
      nodes = t(e,:);
      xy = p(nodes,:);

      Ke = triP1_stiffness(xy);
      fe = triP1_load(xy, f_handle);

      dofs = node2dof(nodes);

      loc = zeros(1,3);
      for a = 1:3
        if dofs(a) > 0
          la = sub(i).glob2loc(dofs(a));
          if la <= 0
            error('assemble_subdomain_matrices_P1: sub(%d) is missing a free DOF that appears in its element set (global dof=%d).', i, dofs(a));
          end
          if la ~= fix(la) || la > nloc
            error('assemble_subdomain_matrices_P1: sub(%d) has invalid local index for a used DOF (global dof=%d).', i, dofs(a));
          end
          loc(a) = la;

          IF(idxF) = la;
          VF(idxF) = fe(a);
          idxF = idxF + 1;
        end
      end

      for a = 1:3
        if dofs(a) == 0, continue; end
        la = loc(a);
        for b = 1:3
          if dofs(b) == 0, continue; end
          lb = loc(b);

          I(idxK) = la;
          J(idxK) = lb;
          V(idxK) = Ke(a,b);
          idxK = idxK + 1;
        end
      end
    end

    I = I(1:idxK-1); J = J(1:idxK-1); V = V(1:idxK-1);
    IF = IF(1:idxF-1); VF = VF(1:idxF-1);

    sub(i).K = sparse(I, J, V, nloc, nloc);
    sub(i).f = sparse(IF, ones(size(IF)), VF, nloc, 1);
  end
end