function ok = run_fetidp_block_integration_test()
%RUN_FETIDP_BLOCK_INTEGRATION_TEST  Runner for key Chapter-4 (FETI-DP) unit tests + block integration.
%
% Runs, in a sensible order:
%   - build_problem_data / setup_fetidp / operator units
%   - solve_fetidp / reconstruct units
%   - test_fetidp_block_integration (this “middle-level” global test)

  ok = true;

  this_file = mfilename('fullpath');
  main_dir  = fileparts(this_file);
  root_dir  = fileparts(main_dir);

  addpath(main_dir);
  if exist('setup_paths','file') == 2
    setup_paths();
  end
  addpath(genpath(fullfile(root_dir, 'tests')));

  clear functions;

  req_tests = { ...
    'test_build_problem_data', ...
    'test_setup_fetidp', ...
    'test_solve_tildeS', ...
    'test_applyA_lambda', ...
    'test_applyM_lambda', ...
    'test_solve_fetidp', ...
    'test_reconstruct_fetidp_solution', ...
    'test_fetidp_block_integration' ...
  };

  for k = 1:numel(req_tests)
    if exist(req_tests{k}, 'file') ~= 2
      error('run_fetidp_block_integration_test:missingTestEntry', ...
            'Expected %s.m to be on path.', req_tests{k});
    end
  end

  try
    for k = 1:numel(req_tests)
      feval(req_tests{k});
      fprintf('[PASS] %s\n', req_tests{k});
    end
  catch err
    ok = false;
    fprintf('[FAIL] FETI-DP test run\n  %s\n', err.message);
    rethrow(err);
  end
end