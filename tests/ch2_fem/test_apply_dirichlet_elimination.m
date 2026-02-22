function test_apply_dirichlet_elimination()
%TEST_APPLY_DIRICHLET_ELIMINATION  Strengthened unit tests for apply_dirichlet_elimination.m

  % Keep tests runnable both under the umbrella runner and standalone. % CHANGED
  if exist('setup_paths', 'file') == 2                                % ADDED
    setup_paths();                                                    % ADDED
  end                                                                 % ADDED

  % ------------------------------------------------------------
  % Test 0: Signature sanity (nargin/nargout)
  % ------------------------------------------------------------
  assert(nargin('apply_dirichlet_elimination') == 3, ...              % ADDED
         sprintf('Expected nargin == 3, got %d.', ...                 % ADDED
                 nargin('apply_dirichlet_elimination')));             % ADDED
  assert(nargout('apply_dirichlet_elimination') == 3, ...             % ADDED
         sprintf('Expected nargout == 3, got %d.', ...                % ADDED
                 nargout('apply_dirichlet_elimination')));            % ADDED

  % ------------------------------------------------------------
  % Test 1: Correct extraction on a small hand-constructed system
  % ------------------------------------------------------------
  K = [ 4 -1  0  0  0;
       -1  4 -1  0  0;
        0 -1  4 -1  0;
        0  0 -1  4 -1;
        0  0  0 -1  4 ];
  F = [10; 20; 30; 40; 50];

  dirichlet_nodes = [2; 5];
  [Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, F, dirichlet_nodes);

  free_expected = [1; 3; 4];
  assert(isequal(free_nodes, free_expected), ...
         'Test 1 failed: free_nodes incorrect.');
  assert(isequal(size(Kff), [3, 3]), ...
         'Test 1 failed: Kff has wrong size.');
  assert(isequal(size(Ff), [3, 1]), ...
         'Test 1 failed: Ff has wrong size.');
  assert(isequal(Kff, K(free_expected, free_expected)), ...
         'Test 1 failed: Kff entries incorrect.');
  assert(isequal(Ff, F(free_expected)), ...
         'Test 1 failed: Ff entries incorrect.');

  % Invariants: free_nodes should be sorted unique column vector in 1..N. % ADDED
  N = size(K,1);                                                      % ADDED
  assert(iscolumn(free_nodes), 'Test 1 failed: free_nodes not column.'); % ADDED
  assert(all(diff(free_nodes) > 0), 'Test 1 failed: free_nodes not strictly increasing.'); % ADDED
  assert(all(free_nodes >= 1 & free_nodes <= N), 'Test 1 failed: free_nodes out of range.'); % ADDED

  % ------------------------------------------------------------
  % Test 2: Robustness to duplicates and unsorted dirichlet_nodes
  % ------------------------------------------------------------
  dirichlet_nodes2 = [5, 2, 2, 5]; % CHANGED: row vector + duplicates
  [Kff2, Ff2, free_nodes2] = apply_dirichlet_elimination(K, F, dirichlet_nodes2);

  assert(isequal(free_nodes2, free_expected), ...
         'Test 2 failed: free_nodes not robust to duplicates/ordering.');
  assert(isequal(Kff2, K(free_expected, free_expected)), ...
         'Test 2 failed: Kff not robust to duplicates/ordering.');
  assert(isequal(Ff2, F(free_expected)), ...
         'Test 2 failed: Ff not robust to duplicates/ordering.');

  % ------------------------------------------------------------
  % Test 3: Empty Dirichlet set => no reduction (identity restriction)
  % ------------------------------------------------------------
  [Kff3, Ff3, free_nodes3] = apply_dirichlet_elimination(K, F, []);    % ADDED
  assert(isequal(free_nodes3, (1:N).'), 'Test 3 failed: free_nodes should be all nodes.'); % ADDED
  assert(isequal(Kff3, K), 'Test 3 failed: Kff should equal K when no Dirichlet nodes.'); % ADDED
  assert(isequal(Ff3, F), 'Test 3 failed: Ff should equal F when no Dirichlet nodes.'); % ADDED

  % ------------------------------------------------------------
  % Test 4: All nodes Dirichlet => must error (no free DOFs remain)
  % ------------------------------------------------------------
  local_assert_throws(@() apply_dirichlet_elimination(K, F, (1:N).'), ... % ADDED
                      'Test 4: all nodes Dirichlet should error');         % ADDED

  % ------------------------------------------------------------
  % Test 5: Wrong-arg-count must error (Octave may throw pre-body)
  % ------------------------------------------------------------
  local_assert_throws(@() apply_dirichlet_elimination(), ...
                      'Test 5a: nargin==0 should error');                 % ADDED
  local_assert_throws(@() apply_dirichlet_elimination(K), ...
                      'Test 5b: nargin==1 should error');                 % ADDED
  local_assert_throws(@() apply_dirichlet_elimination(K, F), ...
                      'Test 5c: nargin==2 should error');                 % ADDED
  local_assert_throws(@() apply_dirichlet_elimination(K, F, [1], 0), ...
                      'Test 5d: nargin==4 should error');                 % ADDED

  % ------------------------------------------------------------
  % Test 6: Shape validation for K and F must error
  % ------------------------------------------------------------
  local_assert_throws(@() apply_dirichlet_elimination(ones(2,3), [1;2], 1), ...
                      'Test 6a: non-square K should error');              % ADDED
  local_assert_throws(@() apply_dirichlet_elimination(eye(3), [1;2], 1), ...
                      'Test 6b: wrong-length F should error');            % ADDED
  local_assert_throws(@() apply_dirichlet_elimination(eye(3), [1,2,3], 1), ...
                      'Test 6c: row-vector F should error');              % ADDED

  % ------------------------------------------------------------
  % Test 7: dirichlet_nodes validity (range/integer/type/finite)
  % ------------------------------------------------------------
  bad_cases = {                                                        % ADDED
    { [2; N+1],  'index > N' },                                         % ADDED
    { [0; 2],    'index < 1' },                                         % ADDED
    { [2; 3.5],  'non-integer' },                                       % ADDED
    { [2; NaN],  'NaN index (must be rejected)' },                      % ADDED
    { [2; Inf],  'Inf index (must be rejected)' },                      % ADDED
    { [2; 1+1i], 'complex index (must be rejected)' },                  % ADDED
    { true,      'logical input (mask/flag, not index list)' }          % ADDED
  };                                                                    % ADDED

  for k = 1:numel(bad_cases)                                            % ADDED
    dn = bad_cases{k}{1};                                               % ADDED
    lbl = bad_cases{k}{2};                                              % ADDED
    local_assert_throws(@() apply_dirichlet_elimination(K, F, dn), ...
                        sprintf('Test 7 case %d: %s', k, lbl));          % ADDED
  end                                                                   % ADDED

  % Char indices can be silently interpreted as ASCII codes.              % ADDED
  % Make N large enough so ASCII values are in-range (e.g., '2' -> 50).   % ADDED
  Kbig = eye(60);                                                       % ADDED
  Fbig = ones(60,1);                                                    % ADDED
  local_assert_throws(@() apply_dirichlet_elimination(Kbig, Fbig, '23'), ... % ADDED
                      'Test 7 char: string indices must error (ASCII hazard)'); % ADDED

  % ------------------------------------------------------------
  % Test 8: Minimal size N=1 sanity
  % ------------------------------------------------------------
  K1 = 2;                                                               % ADDED
  F1 = 3;                                                               % ADDED
  [Kff1, Ff1, free1] = apply_dirichlet_elimination(K1, F1, []);         % ADDED
  assert(isequal(Kff1, K1) && isequal(Ff1, F1) && isequal(free1, 1), ...
         'Test 8 failed: N=1 with empty Dirichlet set should return identity restriction.'); % ADDED
  local_assert_throws(@() apply_dirichlet_elimination(K1, F1, 1), ...
                      'Test 8b: N=1 with node 1 Dirichlet should error'); % ADDED

  fprintf('PASS: test_apply_dirichlet_elimination\n');
end


% ------------------------------------------------------------
% Helper: assert that a function call throws (don’t match message substrings)
% ------------------------------------------------------------
function local_assert_throws(fhandle, label)                            % ADDED
  did_error = false;                                                    % ADDED
  try                                                                   % ADDED
    fhandle();                                                          % ADDED
  catch                                                                 % ADDED
    did_error = true;                                                   % ADDED
  end                                                                   % ADDED
  assert(did_error, sprintf('Expected an error but none was thrown (%s).', label)); % ADDED
end                                                                     % ADDED