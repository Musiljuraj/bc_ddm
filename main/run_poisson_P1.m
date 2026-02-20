% run_poisson_P1.m
% Full FEM solve pipeline for Poisson on unit square with P1 triangles.
% Dirichlet: x=0 and x=1 (homogeneous, handled by elimination).
% Neumann:   y=0 and y=1 (currently homogeneous => no boundary load term).

clear; close all; clc;

% ------------------------------------------------------------
% 0) Paths and graphics toolkit
% ------------------------------------------------------------
% Run Octave from project root, then:
setup_paths();
%addpath(genpath('src'));
%addpath('test');   % optional

% You verified: graphics_toolkit() returns 'qt' => interactive 3D rotation OK.
graphics_toolkit('qt');

% ------------------------------------------------------------
% 1) Problem data
% ------------------------------------------------------------
n = 16;                 % mesh parameter, h = 1/n
f = @(x,y) -1 + 0*x;     % source term f(x,y) (edit as needed)

% ------------------------------------------------------------
% 2) Mesh generation
% ------------------------------------------------------------
[p, t, bnd] = mesh_unit_square_P1(n);
N  = size(p,1);
Ne = size(t,1);

fprintf('Mesh: n=%d, h=%g, N=%d nodes, Ne=%d triangles\n', n, bnd.h, N, Ne);
fprintf('BC: Dirichlet nodes = %d, Neumann edges = %d\n', ...
        numel(bnd.dirichlet_nodes), size(bnd.neumann_edges,1));

% ------------------------------------------------------------
% 3) Assemble global stiffness matrix K and load vector F
% ------------------------------------------------------------
K = assemble_stiffness_P1(p, t);
F = assemble_load_P1(p, t, f);

fprintf('Assembly: nnz(K)=%d, ||F||2=%e, sum(F)=%e\n', nnz(K), norm(F,2), sum(F));

% Optional quick check: symmetry (numerical)
sym_err = norm(K - K.', 'fro') / max(1, norm(K,'fro'));
fprintf('Check: relative symmetry error of K = %e\n', sym_err);

% ------------------------------------------------------------
% 4) Apply homogeneous Dirichlet by elimination (restriction to free DOFs)
% ------------------------------------------------------------
[Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, F, bnd.dirichlet_nodes);

fprintf('Reduced system: %d free DOFs (Kff is %dx%d)\n', ...
        numel(free_nodes), size(Kff,1), size(Kff,2));

% ------------------------------------------------------------
% 5) Solve reduced linear system
% ------------------------------------------------------------
u_free = Kff \ Ff;

% ------------------------------------------------------------
% 6) Reconstruct full solution vector u (Dirichlet nodes are zero)
% ------------------------------------------------------------
u = zeros(N,1);
u(free_nodes) = u_free;

% Diagnostics: reduced residual and basic norms
res2 = norm(Kff*u_free - Ff, 2);
fprintf('Solve: reduced residual ||Kff*u_free - Ff||2 = %e\n', res2);
fprintf('Solution: ||u||2 = %e, max(u) = %e, min(u) = %e\n', ...
        norm(u,2), max(u), min(u));

% ------------------------------------------------------------
% 7) Visualization (interactive 3D “deflected membrane”)
% ------------------------------------------------------------
figure('Name','FEM solution u_h (3D)');
hSurf = trisurf(t, p(:,1), p(:,2), u);
xlabel('x'); ylabel('y'); zlabel('u_h');
title(sprintf('Poisson P1 FEM solution (n=%d, N=%d)', n, N));
view(3);
grid on;
shading interp;
colorbar;

% Enable mouse-driven rotation (Qt supports this)
rotate3d on;

% Optional: adjust aspect ratio for nicer appearance
axis vis3d;

% Optional: show sparsity pattern (useful in thesis debugging/diagnostics)
%figure('Name','Sparsity pattern of K');
%spy(K);
%title(sprintf('Global stiffness matrix K sparsity (n=%d, N=%d)', n, N));
