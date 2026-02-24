function ok = run_ddm_block_integration_test()
%RUN_DDM_BLOCK_INTEGRATION_TEST  Self-contained runner for DDM integration test.

  ok = true;

  this_file = mfilename('fullpath');
  main_dir  = fileparts(this_file);
  root_dir  = fileparts(main_dir);

  test_dir  = fullfile(root_dir, 'tests', 'ch3_ddm');

  addpath(main_dir);
  addpath(test_dir);

  if exist('setup_paths','file') == 2
    setup_paths();
  end

  clear functions;

  if exist('test_ddm_block_integration','file') ~= 2
    error('run_ddm_block_integration_test:missingTestEntry', ...
          'Expected test_ddm_block_integration.m in: %s', test_dir);
  end

  try
    test_ddm_block_integration();
    fprintf('[PASS] test_ddm_block_integration\n');
  catch err
    ok = false;
    fprintf('[FAIL] DDM integration test\n  %s\n', err.message);
    rethrow(err);
  end
end