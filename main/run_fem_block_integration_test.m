function ok = run_fem_block_integration_test()
%RUN_FEM_BLOCK_INTEGRATION_TEST  Self-contained runner (adds paths automatically).
%
% Expected layout:
%   bc_ddm/
%     main/                (this file, run_fem_block_case.m, setup_path.m, ...)
%     tests/ch2_fem/       (test_block_integration.m OR test_fem_block_integration.m)
%
% Usage (from bc_ddm in Octave):
%   addpath(fullfile(pwd,'main'));
%   run_fem_block_integration_test

  ok = true;

  % Locate this file and infer project root (bc_ddm)
  this_file = mfilename('fullpath');
  main_dir  = fileparts(this_file);          % .../bc_ddm/main
  root_dir  = fileparts(main_dir);           % .../bc_ddm

  test_dir  = fullfile(root_dir, 'tests', 'ch2_fem');

  % Add paths needed to find FEM code and integration test
  if exist(main_dir, 'dir') ~= 7
    error('run_fem_block_integration_test:missingMainDir', ...
          'Main directory not found: %s', main_dir);
  end
  addpath(main_dir);

  if exist(test_dir, 'dir') ~= 7
    error('run_fem_block_integration_test:missingTestsDir', ...
          'Tests directory not found: %s', test_dir);
  end
  addpath(test_dir);

  % If you have helpers under tests/ch2_fem subfolders, you can enable this:
  % addpath(genpath(test_dir));

  % Call your existing path setup (if present)
  if exist('setup_paths','file') == 2
    setup_paths();
  elseif exist('setup_path','file') == 2
    setup_path();
  end

  % Ensure fresh function versions after edits
  clear functions;

  % Decide which integration test entry point exists
  has_new_name = (exist('test_fem_block_integration','file') == 2);
  has_old_name = (exist('test_block_integration','file') == 2);

  if ~(has_new_name || has_old_name)
    error('run_fem_block_integration_test:missingTestEntry', ...
          ['Could not find an integration test function on path. ', ...
           'Expected test_fem_block_integration.m or test_block_integration.m in: %s'], test_dir);
  end

  % Run the test
  try
    if has_new_name
      test_fem_block_integration();
      fprintf('[PASS] test_fem_block_integration\n');
    else
      test_block_integration();
      fprintf('[PASS] test_block_integration\n');
    end
  catch err
    ok = false;
    fprintf('[FAIL] integration test\n  %s\n', err.message);
    rethrow(err);
  end
end