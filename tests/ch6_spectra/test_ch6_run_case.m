% ============================================================
% File: tests/ch6_spectra/test_ch6_run_case.m
% ============================================================
function test_ch6_run_case()
%TEST_CH6_RUN_CASE  Smoke/integration test for ch6_run_case on a tiny PDE instance.
%
% Conventions (consistent with Chapter 6 tests):
%   - Robust path setup (works regardless of current working directory)
%   - Concise PASS/FAIL communication
%   - Tiny configuration so the test runs quickly
%
% Purpose:
%   Verify that ch6_run_case:
%     - builds the problem,
%     - runs FETI-DP (first) and BDDC (second),
%     - returns solver stats,
%     - computes spectra (or skips them cleanly based on nmax),
%   without asserting brittle numerical values.

  fprintf('test_ch6_run_case ... ');

  ensure_project_paths_();

  % ----------------------------
  % Tiny case configuration
  % ----------------------------
  cfg = struct();
  cfg.n      = 16;
  cfg.nSubX  = 2;
  cfg.nSubY  = 2;
  cfg.seed   = 1;

  cfg.tol    = 1e-10;
  cfg.maxit  = 200;

  cfg.do_spectra = true;
  cfg.nmax       = 400;    % allow spectra for this tiny case

  cfg.verbose    = false;
  cfg.store_base = false;
  cfg.store_data = false;

  % ----------------------------
  % Run one case
  % ----------------------------
  out = ch6_run_case(cfg);

  % ----------------------------
  % Structural checks (top-level)
  % ----------------------------
  assert(isstruct(out), 'Expected out to be a struct.');
  assert(isfield(out, 'case_id'), 'Expected out.case_id field.');
  assert(isfield(out, 'fetidp'), 'Expected out.fetidp field.');
  assert(isfield(out, 'bddc'),   'Expected out.bddc field.');

  % ----------------------------
  % FETI-DP checks (FIRST)
  % ----------------------------
  assert(isfield(out.fetidp, 'n') && out.fetidp.n >= 1, 'Expected out.fetidp.n >= 1.');
  assert(isfield(out.fetidp, 'stats'), 'Expected out.fetidp.stats field.');
  assert(isstruct(out.fetidp.stats), 'Expected out.fetidp.stats to be a struct.');

  % PCG fields commonly present in your stats
  assert(isfield(out.fetidp.stats, 'iter'),   'Expected out.fetidp.stats.iter field.');
  assert(isfield(out.fetidp.stats, 'relres'), 'Expected out.fetidp.stats.relres field.');

  % Spectra: either computed or skipped cleanly
  assert(isfield(out.fetidp, 'spec_skipped'), 'Expected out.fetidp.spec_skipped field.');
  assert(isfield(out.fetidp, 'spec'),         'Expected out.fetidp.spec field.');
  if ~out.fetidp.spec_skipped
    assert(isfield(out.fetidp.spec, 'chol_ok') && out.fetidp.spec.chol_ok, ...
           'Expected FETI-DP spec.chol_ok == true when not skipped.');
    assert(isfield(out.fetidp.spec, 'eigvals'), 'Expected FETI-DP spec.eigvals.');
    assert(numel(out.fetidp.spec.eigvals) == out.fetidp.n, ...
           'Expected FETI-DP to return exactly n eigenvalues.');
  end

  % ----------------------------
  % BDDC checks (SECOND)
  % ----------------------------
  assert(isfield(out.bddc, 'n') && out.bddc.n >= 1, 'Expected out.bddc.n >= 1.');
  assert(isfield(out.bddc, 'stats'), 'Expected out.bddc.stats field.');
  assert(isstruct(out.bddc.stats), 'Expected out.bddc.stats to be a struct.');

  assert(isfield(out.bddc.stats, 'iter'),   'Expected out.bddc.stats.iter field.');
  assert(isfield(out.bddc.stats, 'relres'), 'Expected out.bddc.stats.relres field.');

  assert(isfield(out.bddc, 'spec_skipped'), 'Expected out.bddc.spec_skipped field.');
  assert(isfield(out.bddc, 'spec'),         'Expected out.bddc.spec field.');
  if ~out.bddc.spec_skipped
    assert(isfield(out.bddc.spec, 'chol_ok') && out.bddc.spec.chol_ok, ...
           'Expected BDDC spec.chol_ok == true when not skipped.');
    assert(isfield(out.bddc.spec, 'eigvals'), 'Expected BDDC spec.eigvals.');
    assert(numel(out.bddc.spec.eigvals) == out.bddc.n, ...
           'Expected BDDC to return exactly n eigenvalues.');
  end

  fprintf('PASS\n');
end


% ============================================================
% Helpers (path setup) — identical approach as other Chapter 6 tests
% ============================================================

function ensure_project_paths_()
%ENSURE_PROJECT_PATHS_  Ensure src/ and tests/ are on path.
  persistent did;
  if ~isempty(did) && did
    return;
  end

  if exist('setup_paths', 'file') == 2
    setup_paths();

    sp = which('setup_paths');
    if ~isempty(sp)
      maindir = fileparts(sp);
      rootdir = fileparts(maindir);
      addpath(genpath(fullfile(rootdir, 'tests')));
      addpath(genpath(fullfile(rootdir, 'main')));
    end

    did = true;
    return;
  end

  % Fallback: locate project root from this test file location
  thisdir  = fileparts(mfilename('fullpath')); % .../tests/ch6_spectra
  testsdir = fileparts(thisdir);               % .../tests
  rootdir  = fileparts(testsdir);              % project root
  maindir  = fullfile(rootdir, 'main');        % .../main

  addpath(maindir);
  setup_paths();
  addpath(genpath(fullfile(rootdir, 'tests')));
  addpath(genpath(fullfile(rootdir, 'main')));

  did = true;
end