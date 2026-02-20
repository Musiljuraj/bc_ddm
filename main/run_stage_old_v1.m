% run_stage.m
% Generic Octave runner for "pipeline up to current file under study".
% Current stage: build_subdomains_structured.m

clear; close all; clc;

% ---------------------------
% CONFIG (edit these)
% ----------------------------
n     = 8;          % mesh resolution used by mesh_unit_square_P1 (e.g., 16 -> 16x16)
nSubX = 2;           % number of subdomains in x
nSubY = 2;           % number of subdomains in y

% RHS (Poisson): -Delta u = f
f = @(x,y) 1 + 0*x;  % constant load (adjust freely)

% Output verbosity
VERBOSE = false;
SHOW_PLOTS = false;

% ----------------------------
% PATHS (minimal)
% Put all needed .m files in the same folder OR add folders here.
% ----------------------------
  addpath(genpath('src'));


% ----------------------------
% CHAPTER 2: mesh + assembly + Dirichlet elimination
% ----------------------------
if VERBOSE, disp('--- Chapter 2: mesh ---'); end
[p, t, bnd] = mesh_unit_square_P1(n);

if VERBOSE
  printf('Mesh: n=%d, nodes=%d, triangles=%d\n', n, size(p,1), size(t,1));
  if isfield(bnd,'dirichlet_nodes')
    printf('Dirichlet nodes: %d\n', numel(bnd.dirichlet_nodes));
  else
    printf('bnd.dirichlet_nodes not found (check mesh generator output).\n');
  end
end

%{
if VERBOSE, disp('--- Chapter 2: global assembly ---'); end
K = assemble_stiffness_P1(p, t);
F = assemble_load_P1(p, t, f);

if VERBOSE
  printf('Global K: %dx%d (nnz=%d)\n', size(K,1), size(K,2), nnz(K));
  printf('Global F: %dx1\n', numel(F));
end


if VERBOSE, disp('--- Chapter 2: Dirichlet elimination ---'); end
% Expected to return reduced matrix/vector and list of free nodes.
% If your function returns different outputs, adjust this line.
[Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, F, bnd.dirichlet_nodes);

if VERBOSE
  printf('Reduced system: nFree=%d\n', numel(free_nodes));
  printf('Reduced Kff: %dx%d (nnz=%d)\n', size(Kff,1), size(Kff,2), nnz(Kff));
end

% Optional: solve reduced system (kept simple; useful for sanity checking)
if VERBOSE, disp('--- (optional) Solve reduced system ---'); end
u_free = Kff \ Ff;
u = zeros(size(p,1),1);
u(free_nodes) = u_free;

if VERBOSE
  printf('Solution stats: ||u||2=%e, min=%e, max=%e\n', norm(u,2), min(u), max(u));
end
%}

% ----------------------------
% CURRENT STAGE: Chapter 3.1 — build_subdomains_structured
% ----------------------------
if VERBOSE, disp('--- Chapter 3.1: build_subdomains_structured ---'); end
% Expected signature (recommended):
[sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);
% If your implementation uses different inputs/outputs, adjust accordingly.
%[sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);


if VERBOSE
  printf('DDM: Nsub=%d (nSubX=%d, nSubY=%d)\n', ddm.Nsub, nSubX, nSubY);

  if isfield(ddm,'nFree')
    printf('DDM: nFree=%d\n', ddm.nFree);
  end

  if isfield(ddm,'free_nodes')
    printf('DDM: free_nodes=%d\n', numel(ddm.free_nodes));
  end

  % Quick per-subdomain summary (keep compact)
  printf('\nPer-subdomain summary:\n');
  printf(' id  (ix,iy)   #elems   #nodes   nloc\n');
  for s = 1:ddm.Nsub
    ix = NaN; iy = NaN;
    if isfield(sub(s),'ix'), ix = sub(s).ix; end
    if isfield(sub(s),'iy'), iy = sub(s).iy; end
    printf('%3d  (%d,%d)   %6d   %6d   %6d\n', ...
      sub(s).id, ix, iy, numel(sub(s).elems), numel(sub(s).nodes), sub(s).nloc);
  end

  % Multiplicity preview: how many subdomains each global free DOF appears in
  if isfield(ddm,'nFree') && ddm.nFree > 0
    mult = zeros(ddm.nFree,1);
    for s = 1:ddm.Nsub
      mult(sub(s).loc2glob) += 1;
    end
    printf('\nSharing preview over global free DOFs:\n');
    printf(' multiplicity==1 : %d\n', sum(mult==1));
    printf(' multiplicity>=2 : %d (candidate interface)\n', sum(mult>=2));
  end
end


% ----------------------------
% CHAPTER 3.2: identify_interface_dofs
% ----------------------------
if VERBOSE, disp('--- Chapter 3.2: identify_interface_dofs ---'); end
% Expected signature (based on your current implementation):
%   [sub, iface] = identify_interface_dofs(sub, ddm)
[sub, iface] = identify_interface_dofs(sub, ddm);


if VERBOSE
  % Common, robust outputs (guarded by isfield)
  if isfield(iface,'nHat')
    printf('Interface (assembled) size nHat = %d\n', iface.nHat);
  end
  if isfield(iface,'glob')
    printf('iface.glob size = %d\n', numel(iface.glob));
  end
  if isfield(iface,'counts')
    c = iface.counts(:);
    printf('Multiplicity summary (iface.counts): min=%d, max=%d\n', min(c), max(c));
    for k = 1:max(c)
      printf('  multiplicity == %d : %d\n', k, sum(c==k));
    end
  else
    % Fallback: recompute multiplicities from sub(s).loc2glob if needed
    if isfield(ddm,'nFree') && ddm.nFree > 0
      mult = zeros(ddm.nFree,1);
      for s = 1:ddm.Nsub
        mult(sub(s).loc2glob) += 1;
      end
      printf('Multiplicity preview (recomputed): max=%d\n', max(mult));
      for k = 1:max(mult)
        printf('  multiplicity == %d : %d\n', k, sum(mult==k));
      end
    end
  end

  % Per-subdomain interior/interface counts if available
  printf('\nPer-subdomain (I/Gamma) summary (after identify_interface_dofs):\n');
  printf(' id    nI    nG    nloc\n');
  for s = 1:ddm.Nsub
    nI = NaN; nG = NaN;
    if isfield(sub(s),'I'),       nI = numel(sub(s).I);       end
    if isfield(sub(s),'Gamma'),   nG = numel(sub(s).Gamma);   end
    if isfield(sub(s),'glob_I'),  nI = numel(sub(s).glob_I);  end
    if isfield(sub(s),'glob_G'),  nG = numel(sub(s).glob_G);  end
    printf('%3d  %5d %5d %6d\n', sub(s).id, nI, nG, sub(s).nloc);
  end
end


% ----------------------------
% CHAPTER 3.2.1: assemble_subdomain_matrices_P1
% ----------------------------
if VERBOSE, disp('--- Chapter 3.2.1: assemble_subdomain_matrices_P1 ---'); end
% Expected signature (from your attached file):
%   sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f_handle)
sub = assemble_subdomain_matrices_P1(p, t, sub, ddm, f);

if VERBOSE
  printf('\nPer-subdomain local matrix/vector summary (after assemble_subdomain_matrices_P1):\n');
  printf(' id   nloc   nnz(K_i)   ||f_i||2\n');
  for s = 1:ddm.Nsub
    nloc = sub(s).nloc;
    if isfield(sub(s),'K'), nnzKi = nnz(sub(s).K); else nnzKi = NaN; end
    if isfield(sub(s),'f'), fi = sub(s).f; else fi = []; end
    if isempty(fi), nfi = NaN; else nfi = norm(full(fi),2); end
    printf('%3d  %5d   %8d   %e\n', sub(s).id, nloc, nnzKi, nfi);
  end

  % Optional quick global consistency hint (not a proof):
  % Sum of local nnz is not comparable to global nnz (different assembly),
  % but you can at least confirm no subdomain has empty K unexpectedly.
  nEmpty = 0;
  for s = 1:ddm.Nsub
    if sub(s).nloc > 0 && (~isfield(sub(s),'K') || nnz(sub(s).K) == 0)
      nEmpty += 1;
    end
  end
  if nEmpty > 0
    printf('NOTE: %d subdomains have nloc>0 but empty K. Check assembly inputs/loc2glob.\n', nEmpty);
  end
end


% ----------------------------
% CHAPTER 3.2.2: extract_subdomain_blocks
% ----------------------------
if VERBOSE, disp('--- Chapter 3.2.2: extract_subdomain_blocks ---'); end
sub = extract_subdomain_blocks(sub);

if VERBOSE
  printf('\nPer-subdomain block summary (after extract_subdomain_blocks):\n');
  printf(' id    nI    nG    size(K_II)    size(K_gg)    ||f_I||2    ||f_g||2\n');
  for s = 1:ddm.Nsub
    nI = NaN; nG = NaN;
    if isfield(sub(s),'nI'), nI = sub(s).nI; end
    if isfield(sub(s),'nG'), nG = sub(s).nG; end

    KII = sub(s).K_II; Kgg = sub(s).K_gg;
    fI  = sub(s).f_I;  fg  = sub(s).f_g;

    printf('%3d  %5d %5d   %4dx%-4d     %4dx%-4d    %e   %e\n', ...
      sub(s).id, nI, nG, size(KII,1), size(KII,2), size(Kgg,1), size(Kgg,2), ...
      norm(full(fI),2), norm(full(fg),2));
  end
end

% ----------------------------
% CHAPTER 3.2.3: setup_local_schur + example application
% ----------------------------
% Local Schur complement on each subdomain:
%   S_i = K_gg - K_gI * K_II^{-1} * K_Ig
% and condensed RHS:
%   rhs_g = f_g - K_gI * K_II^{-1} * f_I
%
% Keep this stage matrix-free by default (assemble_S=false).
opts = struct();
opts.assemble_S = false;
sub = setup_local_schur(sub, opts);

% Concise required summaries (always printed)
printf('\nLocal Schur stage summary:\n');
printf(' sub   nI    nG    K_II_SPD   S_i_dims\n');
for s = 1:ddm.Nsub
  nI = 0; nG = 0;
  if isfield(sub(s),'nI'), nI = sub(s).nI; end
  if isfield(sub(s),'nG'), nG = sub(s).nG; end

  % SPD check: if nI>0, setup_local_schur should have produced a Cholesky factor
  spd = 1;
  if nI > 0
    spd = (isfield(sub(s),'R_II') && ~isempty(sub(s).R_II));
  end
  if ~spd
    error('run_stage: K_II is not SPD (or chol failed) for subdomain %d.', sub(s).id);
  end

  printf('%4d %5d %5d %10d   %4dx%-4d\n', sub(s).id, nI, nG, spd, nG, nG);
end

% Example Schur application (deterministic)
pick = 0;
for s = 1:ddm.Nsub
  if isfield(sub(s),'nG') && sub(s).nG > 0
    pick = s;
    break;
  end
end
if pick == 0
  printf('No subdomain with nG>0 found; skipping example apply_local_schur.\n');
else
  ng = sub(pick).nG;
  x  = (1:ng).';
  y  = apply_local_schur(sub(pick), x);
  printf('Example apply_local_schur: sub %d, x=%dx1 -> y=%dx1, ||y||2=%e\n', ...
         sub(pick).id, numel(x), numel(y), norm(full(y),2));
end

% ----------------------------
% CHAPTER 3.3: build_product_interface  (Product interface space)
% ----------------------------
% Inserted stage per directive: AFTER setup_local_schur (and example apply)
% Do NOT build R, B, or primal logic here.

[sub, prod] = build_product_interface(sub, iface);

% Concise structural summary (always printed)
N       = numel(sub);
ng_hat  = prod.nHat;
ng_prod = prod.nProd;

printf('\n=== Product Interface Stage ===\n');
printf('  #subdomains           : %d\n', N);
printf('  Assembled interface   : ng_hat  = %d\n', ng_hat);
printf('  Product interface     : ng_prod = %d\n', ng_prod);

nG_sum = 0;
for i = 1:N
  % Prefer gamma_glob (created by build_product_interface), fallback to glob_G
  if isfield(sub(i),'gamma_glob')
    nG_i = numel(sub(i).gamma_glob);
  elseif isfield(sub(i),'glob_G')
    nG_i = numel(sub(i).glob_G);
  else
    nG_i = NaN;
  end
  nG_sum = nG_sum + nG_i;
  printf('  sub(%d): nG_i = %d\n', i, nG_i);
end



% ----------------------------
% Optional visuals (off by default)
% ----------------------------
if SHOW_PLOTS

  % --------------------------
  % Figure 1: mesh (2D)
  % --------------------------
  figure(1); clf;
  trimesh(t, p(:,1), p(:,2));
  view(2); axis equal tight;
  title('Mesh (trimesh)');
  xlabel('x'); ylabel('y');
  grid on;
  drawnow;

  % --------------------------
  % Figure 2: subdomain partition (elements colored by subdomain id)
  % --------------------------
  figure(2); clf; hold on;
  z0 = zeros(size(p,1),1);
  for s = 1:ddm.Nsub
    trisurf(t(sub(s).elems,:), p(:,1), p(:,2), z0 + s);
  end
  view(2); axis equal tight;
  title('Subdomain partition (color = subdomain id)');
  xlabel('x'); ylabel('y');
  colorbar;
  shading flat;
  hold off;
  drawnow;

  % --------------------------
  % Figure 3: interface nodes (red) over mesh (black)
  % --------------------------
  figure(3); clf;
  plot(p(:,1), p(:,2), 'k.', 'markersize', 6); hold on;

  % Determine interface global DOFs list robustly
  iface_glob = [];
  if exist('iface','var') && isstruct(iface) && isfield(iface,'glob')
    iface_glob = iface.glob(:);
  elseif exist('iface','var') && isstruct(iface) && isfield(iface,'gHat')
    iface_glob = iface.gHat(:);
  end

  if ~isempty(iface_glob) && isfield(ddm,'dof2node')
    iface_nodes = ddm.dof2node(iface_glob);
    plot(p(iface_nodes,1), p(iface_nodes,2), 'ro', 'markersize', 6, 'linewidth', 1.5);
    legend('all nodes','interface nodes','location','bestoutside');
    title(sprintf('Interface nodes (count = %d)', numel(iface_nodes)));
  else
    legend('all nodes','location','bestoutside');
    title('Interface nodes (not available in iface.*)');
  end

  axis equal tight;
  xlabel('x'); ylabel('y');
  grid on;
  hold off;
  drawnow;

  % --------------------------
  % Figure 4 (optional): solution u (only if computed)
  % --------------------------
  if exist('u','var') && ~isempty(u) && numel(u) == size(p,1)
    figure(4); clf;
    trisurf(t, p(:,1), p(:,2), u);
    title('P1 FEM solution u');
    xlabel('x'); ylabel('y'); zlabel('u');
    view(3); shading interp; colorbar; grid on;
    drawnow;
  end

end

% ----------------------------
% Workspace note (interactive use)
% ----------------------------
% At this point, you can inspect:
%   p,t,bnd,K,F,Kff,Ff,free_nodes,u,
%   sub,ddm,
% and (if computed) mult.
