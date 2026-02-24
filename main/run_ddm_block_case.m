function out = run_ddm_block_case(cfg)
%RUN_DDM_BLOCK_CASE  End-to-end DDM block run with strong diagnostics.
%
% Chain:
%   mesh -> global FEM -> build_subdomains_structured
%   -> identify_interface_dofs -> assemble_subdomain_matrices_P1
%   -> extract_subdomain_blocks -> setup_local_schur (explicit S)
%   -> build_product_interface -> select_primal_dofs -> build_primal_maps
%   -> compare assembled Schur vs global Schur
%
% Required fields:
%   cfg.n, cfg.nSubX, cfg.nSubY
%
% Optional:
%   cfg.f_handle      (default @(x,y) 0)
%   cfg.mutate_mesh   function handle [p,t,bnd] = mutate_mesh(p,t,bnd)
%   cfg.mutate_sub    function handle [sub,ddm] = mutate_sub(sub,ddm)
%
% Output:
%   out struct with all intermediates + diagnostics.

  if nargin ~= 1 || ~isstruct(cfg)
    error('run_ddm_block_case: expected one input struct cfg.');
  end
  req = {'n','nSubX','nSubY'};
  for k = 1:numel(req)
    if ~isfield(cfg, req{k})
      error('run_ddm_block_case: cfg.%s is required.', req{k});
    end
  end

  n     = cfg.n;
  nSubX = cfg.nSubX;
  nSubY = cfg.nSubY;

  f_handle = @(x,y) 0;
  if isfield(cfg,'f_handle') && ~isempty(cfg.f_handle)
    f_handle = cfg.f_handle;
  end
  if ~isa(f_handle,'function_handle')
    error('run_ddm_block_case: cfg.f_handle must be a function handle.');
  end

  % --- Mesh
  [p, t, bnd] = mesh_unit_square_P1(n);

  if isfield(cfg,'mutate_mesh') && ~isempty(cfg.mutate_mesh)
    if ~isa(cfg.mutate_mesh,'function_handle')
      error('run_ddm_block_case: cfg.mutate_mesh must be a function handle.');
    end
    [p, t, bnd] = cfg.mutate_mesh(p, t, bnd);
  end

  % --- Global reference (Dirichlet reduced)
  K  = assemble_stiffness_P1(p, t);
  F  = assemble_load_P1(p, t, f_handle);
  [Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, full(F), bnd.dirichlet_nodes);

  % --- Decomposition
  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);

  if isfield(cfg,'mutate_sub') && ~isempty(cfg.mutate_sub)
    if ~isa(cfg.mutate_sub,'function_handle')
      error('run_ddm_block_case: cfg.mutate_sub must be a function handle.');
    end
    [sub, ddm] = cfg.mutate_sub(sub, ddm);
  end

  % --- Interface bookkeeping
  [sub, iface] = identify_interface_dofs(sub, ddm);

  % --- Local assembly
  sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f_handle);

  % Lift + sum local contributions in ddm DOF ordering (strong consistency check)
  nFree = ddm.nFree;
  Ksum = sparse(nFree, nFree);
  Fsum = zeros(nFree, 1);
  for i = 1:ddm.Nsub
    g = sub(i).loc2glob(:);
    if ~isempty(g)
      Ksum(g,g) = Ksum(g,g) + sub(i).K;
      Fsum(g)   = Fsum(g)   + full(sub(i).f);
    end
  end

  % Permute global Kff/Ff into ddm dof ordering (do NOT assume same ordering)
  nNodes = size(p,1);
  pos_in_free = zeros(nNodes, 1);
  pos_in_free(free_nodes(:)) = (1:numel(free_nodes))';

  idx = pos_in_free(ddm.dof2node(:));
  if any(idx <= 0)
    error('run_ddm_block_case: ddm.dof2node contains nodes not present in free_nodes.');
  end
  Kff_ddm = Kff(idx, idx);
  Ff_ddm  = Ff(idx);

  % --- Block extraction + local Schur (explicit S for testing)
  sub = extract_subdomain_blocks(sub);
  sub = setup_local_schur(sub, struct('assemble_S', true));

  % --- Product interface (also gives gamma_hat per subdomain)
  [sub, prod] = build_product_interface(sub, iface);

  % --- Primal selection + maps
  primal = select_primal_dofs(p, ddm, iface);
  [sub, primal] = build_primal_maps(sub, ddm, iface, prod, primal);

  % --- Global Schur complement computed directly from global Kff (ddm ordering)
  G = iface.glob(:);                       % ordered = hat order
  all = (1:nFree).';
  I = setdiff(all, G);

  K_GG = Kff_ddm(G,G);
  F_G  = Ff_ddm(G);

  if isempty(I)
    S_global = K_GG;
    g_global = F_G;
  else
    K_II = Kff_ddm(I,I);
    K_IG = Kff_ddm(I,G);
    K_GI = Kff_ddm(G,I);
    F_I  = Ff_ddm(I);

    % Use chol for stability (SPD expected)
    R = chol(K_II);
    X = R \ (R' \ K_IG);
    z = R \ (R' \ F_I);

    S_global = K_GG - K_GI * X;
    g_global = F_G  - K_GI * z;
  end

  % --- Assembled-interface Schur from local Schur complements
  nHat = iface.nHat;
  S_hat = sparse(nHat, nHat);
  g_hat = zeros(nHat, 1);
  for i = 1:ddm.Nsub
    h = sub(i).gamma_hat(:);
    if isempty(h), continue; end
    S_hat(h,h) = S_hat(h,h) + sub(i).S;
    g_hat(h)   = g_hat(h)   + sub(i).g;
  end

  % Diagnostics
  diag = struct();

  % (A) local->global assembled matrix match
  diag.rel_Ksum = norm(Ksum - Kff_ddm, 'fro') / max(1, norm(Kff_ddm, 'fro'));
  diag.rel_Fsum = norm(Fsum - Ff_ddm, 2)      / max(1, norm(Ff_ddm, 2));

  % (B) Schur match
  diag.rel_S = norm(S_hat - S_global, 'fro') / max(1, norm(S_global, 'fro'));
  diag.rel_g = norm(g_hat - g_global, 2)     / max(1, norm(g_global, 2));

  % (C) symmetry of assembled Schur
  diag.rel_sym_S = norm(S_hat - S_hat.', 'fro') / max(1, norm(S_hat, 'fro'));

  % Pack output
  out = struct();
  out.p = p; out.t = t; out.bnd = bnd;

  out.K = K; out.F = F;
  out.Kff = Kff; out.Ff = Ff; out.free_nodes = free_nodes;

  out.sub = sub;
  out.ddm = ddm;
  out.iface = iface;
  out.prod = prod;
  out.primal = primal;

  out.Ksum = Ksum;
  out.Fsum = Fsum;
  out.Kff_ddm = Kff_ddm;
  out.Ff_ddm  = Ff_ddm;

  out.S_hat   = S_hat;
  out.g_hat   = g_hat;
  out.S_global = S_global;
  out.g_global = g_global;

  out.diagnostics = diag;
end