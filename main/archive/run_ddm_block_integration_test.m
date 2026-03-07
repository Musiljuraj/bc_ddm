function ok = run_ddm_block_integration_test()
%RUN_DDM_BLOCK_INTEGRATION_TEST  Runner for DDM integration + local-operator unit tests.
%
% Runs, in order:
%   1) test_apply_local_schur
%   2) test_apply_blockdiag_S
%   3) test_ddm_block_integration

  ok = true;

  this_file = mfilename('fullpath');
  main_dir  = fileparts(this_file);
  root_dir  = fileparts(main_dir);

  % Paths
  addpath(main_dir);
  if exist('setup_paths','file') == 2
    setup_paths();
  end
  addpath(genpath(fullfile(root_dir, 'tests')));

  clear functions;

  % Required tests
  req_tests = { 'test_apply_local_schur', 'test_apply_blockdiag_S', 'test_ddm_block_integration' };
  for k = 1:numel(req_tests)
    if exist(req_tests{k}, 'file') ~= 2
      error('run_ddm_block_integration_test:missingTestEntry', ...
            'Expected %s.m to be on path.', req_tests{k});
    end
  end

  try
    test_apply_local_schur();
    fprintf('[PASS] test_apply_local_schur\n');

    test_apply_blockdiag_S();
    fprintf('[PASS] test_apply_blockdiag_S\n');

    test_ddm_block_integration();
    fprintf('[PASS] test_ddm_block_integration\n');
  catch err
    ok = false;
    fprintf('[FAIL] DDM test run\n  %s\n', err.message);
    rethrow(err);
  end
end