% ============================================================
% File: src/common/diagnostics/build_problem_data.m
% ============================================================
function data = build_problem_data(n, nSubX, nSubY, f_handle)
%BUILD_PROBLEM_DATA  Build the minimal “data” struct used by Chapter 4 solvers.
%
% This is intentionally not a general framework. It is a thin, deterministic
% wrapper that reuses Chapter 2–3 building blocks to create:
%   - mesh + global FEM K,F and Dirichlet-reduced Kff,Ff (reference solve)
%   - subdomain decomposition structs
%   - interface bookkeeping (iface/prod), primal selection/maps (primal)
%   - local Schur data per subdomain (explicit S for local splitting)
%
% Inputs:
%   n       : mesh parameter (unit square with (n+1)x(n+1) nodes)
%   nSubX,Y : structured partition counts (must divide n)
%   f_handle: RHS f(x,y) for Poisson; if empty, uses f=1
%
% Output:
%   data struct with fields used downstream by setup_fetidp and drivers.

  if nargin < 4 || isempty(f_handle)
    f_handle = @(x,y) 1.0;
  end

  % --- Chapter 2: mesh + FEM assembly ---
  [p, t, bnd] = mesh_unit_square_P1(n);

  K = assemble_stiffness_P1(p, t);
  F = assemble_load_P1(p, t, f_handle);

  % Dirichlet elimination in GLOBAL dof numbering (P1 scalar dofs at nodes)
  [Kff, Ff, free] = apply_dirichlet_elimination(K, F, bnd.dirichlet_nodes);

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