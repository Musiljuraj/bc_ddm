function ok = run_verified_tests()
%RUN_VERIFIED_TESTS  Run the curated list of verified tests.
%
% Purpose:
%   Runs only tests that have been reviewed and verified as stable.
%
% Behavior:
%   - Sets paths via main/setup_paths.m and tests/ genpath
%   - Clears cached functions (clear functions)
%   - Runs a fixed list of test handles
%   - Prints PASS/FAIL per test
%   - Errors if any test fails
%
% Usage:
%   ok = run_verified_tests();

  ok = true;

  % Locate project root assuming this file is in <root>/main/
  thisfile  = mfilename('fullpath');
  main_dir  = fileparts(thisfile);
  root_dir  = fileparts(main_dir);

  oldpwd = pwd();
  c = onCleanup(@() cd(oldpwd)); %#ok<NASGU>
  cd(root_dir);

  % Ensure paths are available
  addpath(fullfile(root_dir, 'main'));
  setup_paths();
  addpath(genpath(fullfile(root_dir, 'tests')));

  clear functions;

  tests = {
    @test_mesh_unit_square_P1
    @test_triP1_stiffness              % verified leaf FEM element stiffness
    @test_triP1_load                   % verified leaf FEM element load
    @test_apply_dirichlet_elimination  % verified FEM Dirichlet restriction (BC elimination)
    @test_assemble_load_P1             % verified FEM assembly load

    @test_build_subdomains_structured  % verified structured DDM partitioning + mappings
    @test_identify_interface_dofs      % verified DDM interface DOF bookkeeping + maps
    @test_assemble_subdomain_matrices_P1 % verified local subdomain assembly (K,f) vs global Dirichlet-reduced system
    @test_extract_subdomain_blocks     % verified block slicing (K_II,K_Ig,K_gI,K_gg) + RHS splits + input rejection

    @test_setup_local_schur            % verified DDM local Schur setup (R_II,g; nI==0; opts default; input rejection)
    @test_apply_local_schur            % verified local Schur apply (matrix-free) vs explicit Schur + contract checks

    @test_build_product_interface      % verified product interface space bookkeeping (prod2sub/prod2hat/hat2prod) + input rejection
    @test_build_jump_operator_B        % ADDED: verified jump-operator assembly + input rejection
    @test_multiplicity_scaling         % verified multiplicity scaling weights + input rejection
    @test_pcg_wrap                     % verified pcg_wrap contract + solve invariants + input rejection
    @test_rng_deterministic            % verified deterministic RNG seeding + input rejection
    @test_apply_blockdiag_S            % verified block-diagonal Schur apply in product space + input rejection

    @test_select_primal_dofs           % verified primal (corner) interface DOF selection + input contract
    @test_build_primal_maps            % add after verification (maps/splits + hardened contract)

    @test_build_problem_data           % verified end-to-end problem-data builder (tiny DDM config + input rejection)

    @test_setup_fetidp                 % verified FETI-DP setup invariants + input rejection
    @test_applyA_lambda                % verified FETI-DP operator apply wiring + strict input rejection
    @test_applyM_lambda                % verified FETI-DP M_lambda apply (reference equivalence + invariances + input rejection)
    @test_solve_tildeS                 % ADDED: verified solve_tildeS (coarse+delta solve, defaults, linearity, negatives)
    @test_solve_fetidp                 % ADDED: verified solve_fetidp smoke + residual consistency
  };

  failures = {};

  for k = 1:numel(tests)
    name = func2str(tests{k});
    try
      feval(tests{k});
      fprintf('[PASS] %s\n', name);
    catch err
      ok = false;
      fprintf('[FAIL] %s\n  %s\n', name, err.message);
      failures{end+1} = struct('name', name, 'err', err); %#ok<AGROW>
    end
  end

  if ok
    fprintf('\nAll %d verified tests passed.\n', numel(tests));
  else
    fprintf('\n%d/%d verified tests failed.\n', numel(failures), numel(tests));
    error('run_verified_tests:failed', 'Some verified tests failed.');
  end
end