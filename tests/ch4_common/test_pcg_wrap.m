function test_pcg_wrap()
%TEST_PCG_WRAP  Unit tests for pcg_wrap (Chapter 4 common infrastructure).
%
% Run from project root:
%   addpath('main'); setup_paths(); addpath(genpath('tests')); test_pcg_wrap();
%
% Or run from this test's folder directly:
%   test_pcg_wrap();

  ensure_project_paths_();

  % -------------------------
  % Negative tests (must throw)
  % -------------------------
  A = [4, -1; -1, 3];
  applyA = @(x) A * x;
  b = [1; 2];

  bad = {
    @() pcg_wrap()                                                    % wrong arg count
    @() pcg_wrap(applyA, b, 1e-8)                                     % wrong arg count
    @() pcg_wrap(123, b, 1e-8, 10)                                    % applyA not handle
    @() pcg_wrap(applyA, [1 2; 3 4], 1e-8, 10)                        % b not vector
    @() pcg_wrap(applyA, b, -1, 10)                                   % tol negative
    @() pcg_wrap(applyA, b, NaN, 10)                                  % tol not finite
    @() pcg_wrap(applyA, b, 1e-8, 1.5)                                % maxit not integer
    @() pcg_wrap(applyA, b, 1e-8, -2)                                 % maxit negative
    @() pcg_wrap(applyA, b, 1e-8, 10, 7)                              % applyM non-empty non-handle
    @() pcg_wrap(applyA, b, 1e-8, 10, [], ones(3,1))                  % x0 wrong length
    @() pcg_wrap(applyA, b, 1e-8, 10, [], [1 2; 3 4])                 % x0 not vector
  };

  for k = 1:numel(bad)
    assert_throws_(bad{k}, k);
  end

  % ---------------------------------------
  % Positive tests: basic solve, invariants
  % ---------------------------------------
  tol = 1e-12;
  maxit = 50;

  [x, stats] = pcg_wrap(applyA, b, tol, maxit);

  x_ref = A \ b;
  assert_close_(x, x_ref, 1e-10, 1e-12);

  assert(isstruct(stats), 'test_pcg_wrap: stats must be a struct');
  assert(has_field_(stats, 'flag'),   'test_pcg_wrap: stats.flag missing');
  assert(has_field_(stats, 'relres'), 'test_pcg_wrap: stats.relres missing');
  assert(has_field_(stats, 'iter'),   'test_pcg_wrap: stats.iter missing');
  assert(has_field_(stats, 'resvec'), 'test_pcg_wrap: stats.resvec missing');

  assert(isscalar(stats.flag),   'test_pcg_wrap: stats.flag must be scalar');
  assert(isscalar(stats.relres), 'test_pcg_wrap: stats.relres must be scalar');
  assert(isscalar(stats.iter),   'test_pcg_wrap: stats.iter must be scalar');

  % Convergence expectations for a tiny SPD 2x2 system (should be robust).
  assert(stats.flag == 0, 'test_pcg_wrap: expected stats.flag == 0 (converged)');
  assert(stats.relres <= max(tol, 50*eps), 'test_pcg_wrap: relres too large');

  % Residual history consistency
  assert(isvector(stats.resvec), 'test_pcg_wrap: stats.resvec must be a vector');
  assert(numel(stats.resvec) == stats.iter + 1, 'test_pcg_wrap: resvec length must be iter+1');

  % ---------------------------------------
  % Positive: row-vector b should be accepted (wrapper reshapes)
  % ---------------------------------------
  b_row = b.'; % 1x2
  [x2, stats2] = pcg_wrap(applyA, b_row, tol, maxit);
  assert_close_(x2, x_ref, 1e-10, 1e-12);
  assert(stats2.flag == 0, 'test_pcg_wrap: expected convergence for row b case');

  % ---------------------------------------
  % Positive: with a simple Jacobi preconditioner handle
  % ---------------------------------------
  d = diag(A);
  applyM = @(r) r ./ d; % approximate A \ r via diagonal inverse

  [x3, stats3] = pcg_wrap(applyA, b, tol, maxit, applyM);
  assert_close_(x3, x_ref, 1e-10, 1e-12);
  assert(stats3.flag == 0, 'test_pcg_wrap: expected convergence with preconditioner');

  % ---------------------------------------
  % Positive: nonzero initial guess x0
  % ---------------------------------------
  x0 = [10; -3];
  [x4, stats4] = pcg_wrap(applyA, b, tol, maxit, [], x0);
  assert_close_(x4, x_ref, 1e-10, 1e-12);
  assert(stats4.flag == 0, 'test_pcg_wrap: expected convergence with nonzero x0');

  fprintf('PASS: test_pcg_wrap\n');
end

% =========================
% Local helpers (test-only)
% =========================

function ensure_project_paths_()
  persistent did;
  if ~isempty(did) && did
    return;
  end

  if exist('setup_paths', 'file') == 2
    setup_paths();
    sp = which('setup_paths');
    maindir = fileparts(sp);
    rootdir = fileparts(maindir);
    addpath(genpath(fullfile(rootdir, 'tests')));
  else
    thisdir = fileparts(mfilename('fullpath'));
    testsdir = fileparts(thisdir);
    rootdir = fileparts(testsdir);
    maindir = fullfile(rootdir, 'main');
    addpath(maindir);
    setup_paths();
    addpath(genpath(fullfile(rootdir, 'tests')));
  end

  did = true;
end

function assert_throws_(fh, k)
  try
    fh();
  catch
    return;
  end
  error('test_pcg_wrap:expectedError', ...
        'Expected an error but none was thrown (bad-case index %d).', k);
end

function assert_close_(x, y, rtol, atol)
  x = x(:);
  y = y(:);
  if numel(x) ~= numel(y)
    error('test_pcg_wrap:assert_close', 'Size mismatch in assert_close_.');
  end
  dx = norm(x - y);
  ny = norm(y);
  tol = max(atol, rtol * max(1, ny));
  if ~(dx <= tol)
    error('test_pcg_wrap:assert_close', ...
          'Not close: norm(x-y)=%.3e > tol=%.3e', dx, tol);
  end
end

function tf = has_field_(s, name)
  tf = isfield(s, name);
end