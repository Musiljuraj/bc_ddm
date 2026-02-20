function test_pcg_wrap_on_Kff()
%TEST_PCG_WRAP_ON_KFF  Infrastructure test for pcg_wrap using Chapter 2 SPD Kff.
%
% Deterministic, minimal test (Chapter 4.2):
%   - build global FEM stiffness K on the unit square
%   - eliminate Dirichlet DOFs to obtain SPD Kff
%   - run pcg_wrap with applyA = @(x)Kff*x and b = ones
%   - assert successful convergence and basic resvec consistency
%
% Expected environment:
%   The project path setup must make Chapter 2 FEM routines and the Chapter 4.2
%   common tools visible (e.g. via bc_ddm/main/setup_paths.m).

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/fem/elements');
  %addpath('../../src/fem/assembly');
  %addpath('../../src/fem/bc');

  %addpath('../../src/ddm/partition');
  %addpath('../../src/ddm/interface');
  %addpath('../../src/ddm/local');
  %addpath('../../src/ddm/coarse');

  %addpath('../../src/common/coarse');
  %addpath('../../src/common/krylov');
  %addpath('../../src/common/scaling');
  %addpath('../../src/common/utils');

  rng_deterministic(1);

  n = 8;
  [p, t, bnd] = mesh_unit_square_P1(n);

  K = assemble_stiffness_P1(p, t);
  F = ones(size(K,1), 1);
  [Kff, ~, ~] = apply_dirichlet_elimination(K, F, bnd.dirichlet_nodes);

  applyA = @(x) Kff * x;
  b = ones(size(Kff,1), 1);

  tol = 1e-8;
  maxit = 500;
  applyM = [];
  x0 = [];

  [x, stats] = pcg_wrap(applyA, b, tol, maxit, applyM, x0);

  assert(stats.flag == 0, 'pcg_wrap failed to converge (flag=%d).', stats.flag);
  assert(numel(stats.resvec) == stats.iter + 1, 'resvec length mismatch.');
  assert(stats.resvec(end) <= tol * stats.resvec(1), 'residual reduction check failed.');
  assert(isvector(x) && numel(x) == numel(b), 'solution vector size mismatch.');

  fprintf('PASS: test_pcg_wrap_on_Kff\n');
end