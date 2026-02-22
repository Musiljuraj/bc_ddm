function ok = run_verified_tests()
%RUN_VERIFIED_TESTS  Run the curated list of verified tests.
%
% Usage:
%   ok = run_verified_tests();

  ok = true;

  % Locate project root assuming this file is in <root>/main/            % FIXED
  thisfile  = mfilename('fullpath');                                    % CHANGED
  main_dir  = fileparts(thisfile);                                      % ADDED
  root_dir  = fileparts(main_dir);                                      % CHANGED

  oldpwd = pwd();
  c = onCleanup(@() cd(oldpwd));
  cd(root_dir);

  % Ensure paths are available
  addpath(fullfile(root_dir, 'main'));
  setup_paths();
  addpath(genpath(fullfile(root_dir, 'tests')));

  clear functions;

  tests = {
    @test_mesh_unit_square_P1
    @test_triP1_stiffness           % verified leaf FEM element stiffness
    @test_triP1_load                % verified leaf FEM element load      % FIXED
    @test_apply_dirichlet_elimination % ADDED: verified FEM Dirichlet restriction (BC elimination)
    @test_assemble_load_P1          % ADDED: verified FEM assembly load
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