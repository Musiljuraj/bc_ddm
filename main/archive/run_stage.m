% run_stage.m
% Generic Octave runner for "pipeline up to current file under study".
% Current stage: build_product_interface.m
%
% Design rule: main execution is quiet and minimal; all printing goes under
% if VERBOSE ... end, consistent with earlier parts.

clear; close all; clc;

% ---------------------------
% CONFIG (edit these)
% ----------------------------
n     = 8;           % mesh resolution used by mesh_unit_square_P1 (e.g., 16 -> 16x16)
nSubX = 2;           % number of subdomains in x
nSubY = 2;           % number of subdomains in y

% RHS (Poisson): -Delta u = f
f = @(x,y) 1 + 0*x;  % constant load (adjust freely)

% Output verbosity
VERBOSE    = false;
SHOW_PLOTS = false;

% ----------------------------
% PATHS (minimal)
% Put all needed .m files in the same folder OR add folders here.
% ----------------------------
setup_paths();
%addpath(genpath('src'));

% ----------------------------
% CHAPTER 2: mesh + (optional) global FEM solve
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
% --- Optional global assembly/solve (kept disabled by default)
if VERBOSE, disp('--- Chapter 2: global assembly ---'); end
K = assemble_stiffness_P1(p, t);
F = assemble_load_P1(p, t, f);

if VERBOSE
  printf('Global K: %dx%d (nnz=%d)\n', size(K,1), size(K,2), nnz(K));
  printf('Global F: %dx1\n', numel(F));
end

if VERBOSE, disp('--- Chapter 2: Dirichlet elimination ---'); end
[Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, F, bnd.dirichlet_nodes);

if VERBOSE
  printf('Reduced system: nFree=%d\n', numel(free_nodes));
  printf('Reduced Kff: %dx%d (nnz=%d)\n', size(Kff,1), size(Kff,2), nnz(Kff));
end

if VERBOSE, disp('--- (optional) Solve reduced system ---'); end
u_free = Kff \ Ff;
u = zeros(size(p,1),1);
u(free_nodes) = u_free;

if VERBOSE
  printf('Solution stats: ||u||2=%e, min=%e, max=%e\n', norm(u,2), min(u), max(u));
end
%}

% ----------------------------
% CHAPTER 3.1 — build_subdomains_structured
% ----------------------------
if VERBOSE, disp('--- Chapter 3.1: build_subdomains_structured ---'); end
[sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);

if VERBOSE
  printf('DDM: Nsub=%d (nSubX=%d, nSubY=%d)\n', ddm.Nsub, nSubX, nSubY);

  if isfield(ddm,'nFree'),       printf('DDM: nFree=%d\n', ddm.nFree); end
  if isfield(ddm,'free_nodes'),  printf('DDM: free_nodes=%d\n', numel(ddm.free_nodes)); end

  printf('\nPer-subdomain summary:\n');
  printf(' id  (ix,iy)   #elems   #nodes   nloc\n');
  for s = 1:ddm.Nsub
    ix = NaN; iy = NaN;
    if isfield(sub(s),'ix'), ix = sub(s).ix; end
    if isfield(sub(s),'iy'), iy = sub(s).iy; end
    printf('%3d  (%d,%d)   %6d   %6d   %6d\n', ...
      sub(s).id, ix, iy, numel(sub(s).elems), numel(sub(s).nodes), sub(s).nloc);
  end

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
[sub, iface] = identify_interface_dofs(sub, ddm);

if VERBOSE
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
  printf(' id    nI    nG   size(KII)     size(Kgg)    ||fI||2   ||fg||2\n');
  for s = 1:ddm.Nsub
    nI = 0; nG = 0;
    if isfield(sub(s),'glob_I'), nI = numel(sub(s).glob_I); end
    if isfield(sub(s),'glob_G'), nG = numel(sub(s).glob_G); end

    if isfield(sub(s),'K_II'), KII = sub(s).K_II; else KII = sparse(0,0); end
    if isfield(sub(s),'K_gg'), Kgg = sub(s).K_gg; else Kgg = sparse(0,0); end
    if isfield(sub(s),'f_I'),  fI  = sub(s).f_I;  else fI  = zeros(0,1); end
    if isfield(sub(s),'f_g'),  fg  = sub(s).f_g;  else fg  = zeros(0,1); end

    printf('%3d  %5d %5d   %4dx%-4d     %4dx%-4d    %e   %e\n', ...
      sub(s).id, nI, nG, size(KII,1), size(KII,2), size(Kgg,1), size(Kgg,2), ...
      norm(full(fI),2), norm(full(fg),2));
  end
end

% ----------------------------
% CHAPTER 3.2.3: setup_local_schur (+ optional example application under VERBOSE)
% ----------------------------
% Local Schur complement on each subdomain:
%   S_i = K_gg - K_gI * K_II^{-1} * K_Ig
% and condensed RHS:
%   rhs_g = f_g - K_gI * K_II^{-1} * f_I
%
% Keep this stage matrix-free by default (assemble_S=false).
sub = setup_local_schur(sub, struct('assemble_S', false));

% Quiet always-on consistency check (no printing, minimal locals).
% If nI>0 then setup_local_schur is expected to have produced a Cholesky factor.
for s = 1:ddm.Nsub
  nI = 0;
  if isfield(sub(s),'nI'), nI = sub(s).nI; end
  if nI > 0
    if ~(isfield(sub(s),'R_II') && ~isempty(sub(s).R_II))
      error('run_stage: K_II is not SPD (or chol failed) for subdomain %d.', sub(s).id);
    end
  end
end

if VERBOSE
  printf('\nLocal Schur stage summary:\n');
  printf(' sub   nI    nG    K_II_SPD   S_i_dims\n');
  for s = 1:ddm.Nsub
    nI = 0; nG = 0;
    if isfield(sub(s),'nI'), nI = sub(s).nI; end
    if isfield(sub(s),'nG'), nG = sub(s).nG; end
    spd = 1;
    if nI > 0
      spd = (isfield(sub(s),'R_II') && ~isempty(sub(s).R_II));
    end
    printf('%4d %5d %5d %10d   %4dx%-4d\n', sub(s).id, nI, nG, spd, nG, nG);
  end

  % Deterministic example Schur application (kept under VERBOSE)
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
end

% ----------------------------
% CHAPTER 3.3: build_product_interface  (Product interface space)
% ----------------------------
[sub, prod] = build_product_interface(sub, iface);

% Quiet always-on light checks (no printing, minimal locals)
nG_sum = 0;
for i = 1:numel(sub)
  if isfield(sub(i),'gamma_glob')
    nG_sum += numel(sub(i).gamma_glob);
  elseif isfield(sub(i),'glob_G')
    nG_sum += numel(sub(i).glob_G);
  end
end

if prod.nProd ~= nG_sum
  error('run_stage: product interface size mismatch: nProd=%d, sum_i nG_i=%d.', prod.nProd, nG_sum);
end

if ~all(prod.prod2hat >= 1 & prod.prod2hat <= prod.nHat)
  error('run_stage: prod2hat has entries outside [1, nHat].');
end

if VERBOSE
  printf('\n=== Product Interface Stage ===\n');
  printf('  #subdomains           : %d\n', numel(sub));
  printf('  Assembled interface   : ng_hat  = %d\n', prod.nHat);
  printf('  Product interface     : ng_prod = %d\n', prod.nProd);

  nG_sum_v = 0;
  for i = 1:numel(sub)
    if isfield(sub(i),'gamma_glob')
      nG_i = numel(sub(i).gamma_glob);
    elseif isfield(sub(i),'glob_G')
      nG_i = numel(sub(i).glob_G);
    else
      nG_i = NaN;
    end
    nG_sum_v += nG_i;
    printf('  sub(%d): nG_i = %d\n', i, nG_i);
  end

  printf('  [OK] ng_prod == sum_i nG_i\n');
  printf('  [OK] prod2hat indices within [1, ng_hat]\n');

  % Deterministic sanity test (no randomness)
  u_hat = (1:prod.nHat).';
  u_prod = u_hat(prod.prod2hat);

  % Explicit gathered reference (same definition), for a norm check
  u_ref = zeros(prod.nProd, 1);
  for j = 1:prod.nProd
    u_ref(j) = u_hat(prod.prod2hat(j));
  end
  printf('  sanity: norm(u_prod - gathered_reference) = %.6e\n', norm(u_prod - u_ref));
  printf('=== End Product Interface Stage ===\n\n');
end


% ----------------------------
% CHAPTER 3.4: build_assembly_operator_R  (Assembly operator R: hat -> product)
% ----------------------------
R = build_assembly_operator_R(prod);

% Quiet always-on structural checks (fail fast)
if size(R,1) ~= prod.nProd
  error('run_stage: R has wrong number of rows: size(R,1)=%d, expected nProd=%d.', size(R,1), prod.nProd);
end
if size(R,2) ~= prod.nHat
  error('run_stage: R has wrong number of cols: size(R,2)=%d, expected nHat=%d.', size(R,2), prod.nHat);
end

% Canonical unweighted gather/distribution: exactly one nonzero per product DOF row
if nnz(R) ~= prod.nProd
  error('run_stage: R has unexpected nnz(R)=%d (expected nProd=%d for canonical 1-per-row R).', nnz(R), prod.nProd);
end

% Deterministic mapping sanity check (no workspace-cluttering intermediates)
if norm(R * (1:prod.nHat).' - ((1:prod.nHat).'(prod.prod2hat)), inf) ~= 0
  error('run_stage: R mapping check failed: R*u_hat != u_hat(prod2hat) for u_hat=(1:nHat).');
end

if VERBOSE
  printf('\n=== Assembly Operator R Stage ===\n');
  printf('  size(R) : %d x %d\n', size(R,1), size(R,2));
  printf('  nnz(R)  : %d\n', nnz(R));
  printf('  mapping : OK (R*(1:nHat)'' matches direct gather via prod2hat)\n');
  printf('=== End Assembly Operator R Stage ===\n\n');
end

% ----------------------------
% CHAPTER 3.5: build_jump_operator_B  (Jump/constraint operator B in product space)
% ----------------------------
B = build_jump_operator_B(prod);

% Quiet always-on structural checks (fail fast)
if size(B,2) ~= prod.nProd
  error('run_stage: B has wrong number of cols: size(B,2)=%d, expected nProd=%d.', size(B,2), prod.nProd);
end
if size(R,1) ~= prod.nProd
  error('run_stage: R has wrong number of rows: size(R,1)=%d, expected nProd=%d.', size(R,1), prod.nProd);
end

% Algebraic identity check: B*R == 0  (Range(R) subset Ker(B))
M = B * R;

if size(M,2) ~= prod.nHat
  error('run_stage: B*R has wrong number of cols: size(B*R,2)=%d, expected nHat=%d.', size(M,2), prod.nHat);
end

tol_BR = 1e-14;
if norm(M, inf) > tol_BR
  error('run_stage: identity check failed: norm(B*R,inf)=%.3e > %.3e.', norm(M, inf), tol_BR);
end

% Optional light sanity check: each constraint row should sum to zero
% (balanced +/- entries). Skip if B has zero rows.
if size(B,1) > 0
  if norm(sum(B,2), inf) > tol_BR
    error('run_stage: sanity check failed: some rows of B do not sum to zero (inf-norm %.3e).', norm(sum(B,2), inf));
  end
end

if VERBOSE
  printf('\n=== Jump Operator B Stage ===\n');
  printf('  size(B) : %d x %d\n', size(B,1), size(B,2));
  printf('  nConstr : %d\n', size(B,1));
  printf('  check   : OK (norm(B*R,inf) <= %.1e)\n', tol_BR);
  if size(B,1) > 0
    printf('  rowsum  : OK (each row sums to zero)\n');
  end
  printf('=== End Jump Operator B Stage ===\n\n');
end

% ----------------------------
% CHAPTER 3.6: apply_blockdiag_S  (Block-diagonal Schur operator in product space)
% ----------------------------
% Deterministic product-space test vector (no randomness)
w = (1:prod.nProd).';

% Apply block-diagonal Schur operator in product space
% NOTE: function signature is apply_blockdiag_S(sub, w)
y = apply_blockdiag_S(sub, w);

% Minimal fail-fast checks
if numel(y) ~= prod.nProd || size(y,2) ~= 1
  error('run_stage: apply_blockdiag_S: y must be a column vector of length nProd.');
end
if any(~isfinite(y))
  error('run_stage: apply_blockdiag_S: y contains NaN/Inf.');
end

% Strong local consistency check (cheap, blockwise)
tol_S = 1e-12;
for i = 1:numel(sub)
  idx = sub(i).prod_idx(:);
  if isempty(idx), continue; end
  if norm(y(idx) - apply_local_schur(sub(i), w(idx)), inf) > tol_S
    error('run_stage: apply_blockdiag_S mismatch on sub(%d) (inf-norm > %.1e).', i, tol_S);
  end
end

if VERBOSE
  printf('\n=== Block-Diagonal Schur Stage ===\n');
  printf('  check   : OK (blockwise match to apply_local_schur, tol=%.1e)\n', tol_S);
  printf('=== End Block-Diagonal Schur Stage ===\n\n');
end


% ----------------------------
% CHAPTER 3.7: select_primal_dofs  (Primal/coarse DOF selection on the interface)
% ----------------------------
primal = select_primal_dofs(p, ddm, iface);

% Minimal structural checks (fail fast)
if ~isstruct(primal) || ~isfield(primal,'glob_c') || ~isfield(primal,'hat_c') || ~isfield(primal,'nC')
  error('run_stage: select_primal_dofs returned invalid output (missing glob_c/hat_c/nC).');
end

% No duplicates in primal set
if numel(unique(primal.glob_c(:))) ~= numel(primal.glob_c(:))
  error('run_stage: primal.glob_c contains duplicates.');
end

% Primal must be subset of interface DOFs (global numbering)
if ~all(ismember(primal.glob_c(:), iface.glob(:)))
  error('run_stage: primal.glob_c contains non-interface DOFs.');
end

% Hat indices must be in valid assembled interface range
if ~all(primal.hat_c(:) >= 1 & primal.hat_c(:) <= prod.nHat)
  error('run_stage: primal.hat_c has entries outside [1, prod.nHat].');
end

if VERBOSE
  printf('\n=== Primal Selection Stage ===\n');
  printf('  nC (primal DOFs): %d\n', primal.nC);
  if primal.nC > 0
    show = min(10, primal.nC);
    printf('  glob_c (first %d): ', show);
    fprintf('%d ', primal.glob_c(1:show));
    printf('\n');
  end
  printf('=== End Primal Selection Stage ===\n\n');
end

% ----------------------------
% CHAPTER 3.7: select_primal_dofs  (Primal/coarse DOF selection on the interface)
% ----------------------------
primal = select_primal_dofs(p, ddm, iface);

% Minimal checks (keep it light)
if ~isstruct(primal) || ~isfield(primal,'glob_c')
  error('run_stage: select_primal_dofs returned invalid output (missing glob_c).');
end

% ----------------------------
% CHAPTER 3.8: build_primal_maps  (Primal/delta maps & bookkeeping)
% ----------------------------
[sub, primal] = build_primal_maps(sub, ddm, iface, prod, primal);

% Minimal checks: key maps exist + sizes match index spaces
if ~isfield(primal,'glob2c') || numel(primal.glob2c) ~= ddm.nFree
  error('run_stage: build_primal_maps produced invalid glob2c.');
end
if ~isfield(primal,'hat2c') || numel(primal.hat2c) ~= iface.nHat
  error('run_stage: build_primal_maps produced invalid hat2c.');
end
if ~isfield(primal,'prod_is_c') || numel(primal.prod_is_c) ~= prod.nProd
  error('run_stage: build_primal_maps produced invalid prod_is_c.');
end
if ~isfield(primal,'prod2delta') || numel(primal.prod2delta) ~= prod.nProd
  error('run_stage: build_primal_maps produced invalid prod2delta.');
end

if VERBOSE
  printf('\n=== Primal Maps Stage ===\n');
  if isfield(primal,'nC'),        printf('  nC        : %d\n', primal.nC); end
  if isfield(primal,'nDeltaProd'),printf('  nDeltaProd: %d\n', primal.nDeltaProd); end
  printf('=== End Primal Maps Stage ===\n\n');
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
  title('Mesh (P1 triangles)');

  % --------------------------
  % Figure 2: subdomain partition (by element id)
  % --------------------------
  figure(2); clf; hold on;
  % Plot elements colored by subdomain assignment if available
  if isfield(sub(1),'elems')
    for s = 1:numel(sub)
      elems = sub(s).elems(:);
      for ee = elems.'
        tri = t(ee,:);
        patch(p(tri,1), p(tri,2), s, 'EdgeColor', 'k');
      end
    end
    axis equal tight;
    title('Subdomain partition (element coloring by subdomain index)');
    colorbar;
  else
    title('Subdomain partition: sub(s).elems not found');
  end

end
