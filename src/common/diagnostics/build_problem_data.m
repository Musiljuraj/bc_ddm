% ============================================================
% File: src/common/diagnostics/build_problem_data.m
% ============================================================
function data = build_problem_data(n, nSubX, nSubY, f_handle)
%BUILD_PROBLEM_DATA Build the common data structure for Octave experiments.
% Thesis link: Chapters 3–5 (FEM data, substructuring data, solver data).
% This routine gathers the model, discretization, and decomposition objects
% required by the sequential FETI-DP and BDDC workflows.
%
% Inputs:
%   n       : mesh parameter (unit square with (n+1)x(n+1) nodes)
%   nSubX,Y : structured partition counts (must divide n)
%   f_handle: RHS f(x,y) for Poisson; if empty, uses f=1
%
% Output:
%   data struct with fields used downstream by setup_fetidp and drivers.

  narginchk(3, 4);                                

  validate_pos_int_scalar_(n, 'n');               
  validate_pos_int_scalar_(nSubX, 'nSubX');         
  validate_pos_int_scalar_(nSubY, 'nSubY');          
  if (nSubX * nSubY) < 2                             
    error('build_problem_data:singleSubdomain', ...
          'build_problem_data requires at least 2 subdomains (nSubX*nSubY >= 2).');
  end

  if mod(n, nSubX) ~= 0 || mod(n, nSubY) ~= 0        
    error('build_problem_data:invalidPartition', ...
          'Structured partition requires nSubX and nSubY to divide n.');
  end

  if nargin < 4 || isempty(f_handle)
    f_handle = @(x,y) 1.0;
  else
    if ~isa(f_handle, 'function_handle')             
      error('build_problem_data:invalidRHS', ...
            'f_handle must be a function handle f(x,y) or empty.');
    end
  end

  % --- Chapter 2: mesh + FEM assembly ---
  [p, t, bnd] = mesh_unit_square_P1(n);

  % ADDED: boundary struct validation (fail fast, clearer error than deep stack)
  if ~isstruct(bnd) || ~isfield(bnd, 'dirichlet_nodes')
    error('build_problem_data:missingBoundaryField', ...
          'mesh_unit_square_P1 must return bnd.dirichlet_nodes.');
  end

  % ADDED: validate dirichlet_nodes indices now, before elimination
  nNodes = size(p, 1);
  dirn = bnd.dirichlet_nodes;
  validate_index_vector_(dirn, nNodes, 'bnd.dirichlet_nodes');
  bnd.dirichlet_nodes = unique(dirn(:));             

  K = assemble_stiffness_P1(p, t);
  F = assemble_load_P1(p, t, f_handle);

  % Dirichlet elimination in GLOBAL dof numbering (P1 scalar dofs at nodes)
  [Kff, Ff, free] = apply_dirichlet_elimination(K, F, bnd.dirichlet_nodes);

  if isempty(free)                                   
    error('build_problem_data:emptyFreeDofs', ...
          'Dirichlet elimination produced no free DOFs; increase n.');
  end

  data = struct();
  data.p = p;
  data.t = t;
  data.bnd = bnd;                                    

  data.K = K;
  data.F = F;
  data.Kff = Kff;
  data.Ff = Ff;
  data.free = free;

  % --- Chapter 3: structured decomposition + interface ---
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);

  % Identify interface DOFs in the FREE-DOF numbering used after Dirichlet elimination.
  [sub, iface] = identify_interface_dofs(sub, ddm);

  % ADDED: enforce downstream contract explicitly (select_primal_dofs requires nonempty iface.glob)
  if ~isstruct(iface) || ~isfield(iface, 'glob') || isempty(iface.glob)
    error('build_problem_data:emptyInterface', ...
          'Decomposition produced an empty interface (iface.glob). Require >=2 subdomains and a nontrivial mesh.');
  end

  % Assemble subdomain matrices in FREE-DOF numbering.
  sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f_handle);

  % Extract local blocks (I/G split, etc.)
  sub = extract_subdomain_blocks(sub);

  % Setup local Schur objects; request explicit local Schur S for later splitting.
  opts = struct();
  opts.assemble_S = true;
  sub = setup_local_schur(sub, opts);

  % Product interface bookkeeping
  [sub, prod] = build_product_interface(sub, iface);

  % Primal selection (corners-only) + primal maps/splits
  primal = select_primal_dofs(p, ddm, iface);
  [sub, primal] = build_primal_maps(sub, ddm, iface, prod, primal);

  % Attach
  data.ddm = ddm;
  data.iface = iface;
  data.sub = sub;
  data.prod = prod;
  data.primal = primal;
end

% ========================= local validation helpers =========================
% NOTE: These helpers are local to this file only (no project-wide side effects).

function validate_pos_int_scalar_(x, name)            
  if ~(isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x))
    error('build_problem_data:invalidInput', ...
          '%s must be a real, finite numeric scalar.', name);
  end
  if x <= 0
    error('build_problem_data:invalidInput', ...
          '%s must be positive.', name);
  end
  if x ~= round(x)
    error('build_problem_data:invalidInput', ...
          '%s must be integer-valued.', name);
  end
end

function validate_index_vector_(idx, nMax, name)      
  if isempty(idx)
    error('build_problem_data:invalidIndex', ...
          '%s must be a nonempty vector of indices.', name);
  end
  if ~(isnumeric(idx) && isreal(idx))
    error('build_problem_data:invalidIndex', ...
          '%s must be a real numeric vector.', name);
  end
  idx = idx(:);
  if any(~isfinite(idx))
    error('build_problem_data:invalidIndex', ...
          '%s must contain only finite values.', name);
  end
  if any(idx ~= round(idx))
    error('build_problem_data:invalidIndex', ...
          '%s must be integer-valued.', name);
  end
  if any(idx < 1) || any(idx > nMax)
    error('build_problem_data:invalidIndex', ...
          '%s indices must be in the range 1..%d.', name, nMax);
  end
end