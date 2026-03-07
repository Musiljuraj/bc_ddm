function ok = run_bddc_block_integration_test()
%RUN_BDDC_BLOCK_INTEGRATION_TEST  Runner for BDDC unit tests + block integration.

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
    'test_setup_bddc', ...
    'test_applyA_hat', ...
    'test_applyM_bddc', ...
    'test_solve_bddc', ...
    'test_reconstruct_bddc_solution', ...
    'test_bddc_block_integration' ...
  };

  for k = 1:numel(req_tests)
    if exist(req_tests{k}, 'file') ~= 2
      error('run_bddc_block_integration_test:missingTestEntry', ...
            'Expected %s.m to be on path.', req_tests{k});
    end
  end

  try
    test_setup_bddc();                fprintf('[PASS] test_setup_bddc\n');
    test_applyA_hat();                fprintf('[PASS] test_applyA_hat\n');
    test_applyM_bddc();               fprintf('[PASS] test_applyM_bddc\n');
    test_solve_bddc();                fprintf('[PASS] test_solve_bddc\n');
    test_reconstruct_bddc_solution(); fprintf('[PASS] test_reconstruct_bddc_solution\n');
    test_bddc_block_integration();    fprintf('[PASS] test_bddc_block_integration\n');
  catch err
    ok = false;
    fprintf('[FAIL] BDDC test run\n  %s\n', err.message);
    rethrow(err);
  end
end