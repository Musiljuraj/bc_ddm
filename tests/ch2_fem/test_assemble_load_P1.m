function test_assemble_load_P1()
%TEST_ASSEMBLE_LOAD_P1  Strengthened unit tests for assemble_load_P1.m

  setup_paths();  % (safe even if caller already ran it)  % CHANGED

  % -----------------------------
  % Test 0: wrong argument counts
  % -----------------------------
  p1 = [0 0; 1 0; 0 1];
  t1 = [1 2 3];
  f1 = @(x,y) 1.0;

  assert_throws(@() assemble_load_P1(),                'nargin=0');   % ADDED
  assert_throws(@() assemble_load_P1(p1, t1),          'nargin=2');   % ADDED
  assert_throws(@() assemble_load_P1(p1, t1, f1, 7),   'nargin=4');   % ADDED

  % ----------------------------------------------------
  % Test 1: single triangle, constant forcing (exact)
  % ----------------------------------------------------
  c = 2.5;
  f = @(x,y) c;

  F = assemble_load_P1(p1, t1, f);

  assert(all(size(F) == [3, 1]), ...
         sprintf('Test 1 failed: size(F)=[%d,%d], expected [3,1].', size(F,1), size(F,2))); % CHANGED
  assert(issparse(F), 'Test 1 failed: F should be sparse.');                                % ADDED

  A = tri_area(p1);                     % area = 0.5
  F_expected = (c * A / 3) * ones(3,1); % exact for constant f with P1 basis + centroid rule
  assert_close(full(F), F_expected, 'Test 1 failed: constant load mismatch on single triangle.'); % ADDED

  % ----------------------------------------------------
  % Test 2: orientation / local ordering invariance
  % For constant f, should be identical.
  % ----------------------------------------------------
  t1_flip = [1 3 2];
  F_flip = assemble_load_P1(p1, t1_flip, f);
  assert_close(full(F_flip), full(F), ...
               'Test 2 failed: flipping triangle vertex order changed F for constant f.'); % ADDED

  % ----------------------------------------------------
  % Test 3: two-triangle unit square, constant forcing (exact)
  % ----------------------------------------------------
  p2 = [0 0; 1 0; 0 1; 1 1];
  t2 = [1 2 3; 2 4 3];
  c2 = -3.0;
  f2 = @(x,y) c2;

  F2 = assemble_load_P1(p2, t2, f2);

  assert(all(size(F2) == [4, 1]), 'Test 3 failed: wrong output size for two-triangle mesh.'); % ADDED
  assert(issparse(F2), 'Test 3 failed: F2 should be sparse.');                                % ADDED

  a = 0.5;               % each triangle area
  v = c2 * a / 3;
  F2_expected = [v; (v+v); (v+v); v];  % nodes 2 and 3 shared
  assert_close(full(F2), F2_expected, 'Test 3 failed: constant load mismatch on two triangles.'); % ADDED

  % Global integral check: sum(F) = c2 * area(domain) = c2 * 1
  assert_close(sum(full(F2)), c2, 'Test 3 failed: sum(F2) should equal c2 for unit area.'); % ADDED

  % ----------------------------------------------------
  % Test 4: empty mesh (Ne=0) => zero vector
  % ----------------------------------------------------
  t_empty = zeros(0,3);
  F0 = assemble_load_P1(p2, t_empty, f2);
  assert(all(size(F0) == [size(p2,1), 1]), 'Test 4 failed: empty mesh returned wrong size.'); % ADDED
  assert(issparse(F0), 'Test 4 failed: empty mesh should return sparse vector.');            % ADDED
  assert(nnz(F0) == 0, 'Test 4 failed: empty mesh should produce zero vector.');             % ADDED

  % ----------------------------------------------------
  % Test 5: renumbering invariance (global permutation)
  % ----------------------------------------------------
  perm = [3 1 4 2];
  invperm = zeros(1, numel(perm));
  invperm(perm) = 1:numel(perm);

  p2p = p2(perm, :);
  t2p = invperm(t2);

  F2p = assemble_load_P1(p2p, t2p, f2);
  assert_close(full(F2p), full(F2(perm)), 'Test 5 failed: permutation invariance violated.'); % ADDED

  % ----------------------------------------------------
  % Test 6: linearity in f for constant sources
  % ----------------------------------------------------
  ca = 1.25;
  cb = -0.75;
  fa = @(x,y) ca;
  fb = @(x,y) cb;
  fab = @(x,y) (ca + cb);

  Fa  = assemble_load_P1(p2, t2, fa);
  Fb  = assemble_load_P1(p2, t2, fb);
  Fab = assemble_load_P1(p2, t2, fab);

  assert_close(full(Fab), full(Fa) + full(Fb), 'Test 6 failed: linearity in f violated.'); % ADDED

  % ----------------------------------------------------
  % Test 7: input validation must throw (no message matching)
  % ----------------------------------------------------
  % 7a) p wrong shape
  assert_throws(@() assemble_load_P1([0 0 0], t1, f1), 'p is N x 3'); % ADDED

  % 7b) t wrong shape
  assert_throws(@() assemble_load_P1(p1, [1 2], f1), 't is Ne x 2');  % ADDED

  % 7c) f_handle wrong type
  assert_throws(@() assemble_load_P1(p1, t1, 7), 'f_handle numeric'); % ADDED

  % 7d) t has out-of-range / non-integer indices
  assert_throws(@() assemble_load_P1(p1, [0 2 3], f1),   't contains 0');          % ADDED
  assert_throws(@() assemble_load_P1(p1, [1 2 4], f1),   't contains >N');         % ADDED
  assert_throws(@() assemble_load_P1(p1, [1 2 2.2], f1), 't contains non-integer');% ADDED

  % 7e) NaN/Inf in connectivity should throw (prefer early validation)
  assert_throws(@() assemble_load_P1(p1, [1 NaN 3], f1), 't contains NaN'); % ADDED
  assert_throws(@() assemble_load_P1(p1, [1 Inf 3], f1), 't contains Inf'); % ADDED

  % 7f) p non-numeric / logical / non-finite / complex should throw
  p_char = ['a' 'b'; 'c' 'd'; 'e' 'f'];
  assert_throws(@() assemble_load_P1(p_char, t1, f1), 'p is char'); % ADDED

  p_logical = logical([0 0; 1 0; 0 1]);
  assert_throws(@() assemble_load_P1(p_logical, t1, f1), 'p is logical'); % ADDED

  p_nan = p1; p_nan(1,1) = NaN;
  assert_throws(@() assemble_load_P1(p_nan, t1, f1), 'p contains NaN'); % ADDED

  p_inf = p1; p_inf(2,2) = Inf;
  assert_throws(@() assemble_load_P1(p_inf, t1, f1), 'p contains Inf'); % ADDED

  p_cplx = p1; p_cplx(3,1) = 1 + 1i;
  assert_throws(@() assemble_load_P1(p_cplx, t1, f1), 'p contains complex'); % ADDED

  % 7g) IMPORTANT silent-acceptance cases for t:
  %     - char matrices (ASCII indices)
  %     - logical ones (all indices = 1)
  % These MUST throw; if they do not, that's a real weakness in assemble_load_P1.
  Nbig = 200;
  pbig = [(1:Nbig).' , mod((1:Nbig).',2)];
  t_char = 'abc';  % ASCII indices [97 98 99]
  assert_throws(@() assemble_load_P1(pbig, t_char, f1), 't is char (ASCII indices)'); % ADDED

  t_log = logical([1 1 1]); % would silently assemble into node 1 otherwise
  assert_throws(@() assemble_load_P1(p1, t_log, f1), 't is logical'); % ADDED

  % ----------------------------------------------------
  % Test 8: pipeline sanity (mesh_unit_square_P1) like your original
  % ----------------------------------------------------
  n = 4;
  [p, t, bnd] = mesh_unit_square_P1(n); %#ok<ASGLU>
  c3 = 2.5;
  f3 = @(x,y) c3;
  F3 = assemble_load_P1(p, t, f3);

  assert(all(size(F3) == [size(p,1), 1]), 'Test 8 failed: wrong size on unit square mesh.'); % ADDED
  assert_close(sum(full(F3)), c3, 'Test 8 failed: sum(F) should equal c on unit square mesh.'); % ADDED

  fprintf('PASS: test_assemble_load_P1\n');
end

% =========================
% Helper utilities (local)
% =========================

function A = tri_area(xy) % ADDED
  x1 = xy(1,1); y1 = xy(1,2);
  x2 = xy(2,1); y2 = xy(2,2);
  x3 = xy(3,1); y3 = xy(3,2);
  A = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
end

function assert_close(a, b, label) % ADDED
  a = a(:); b = b(:);
  scale = max(1, max(abs(b)));
  tol = 200 * eps * scale;
  err = norm(a - b, inf);
  assert(err <= tol, sprintf('%s (err=%g, tol=%g)', label, err, tol));
end

function assert_throws(fun, label) % ADDED (Rule 7.3 diagnostic)
  threw = false;
  try
    fun();
  catch
    threw = true;
  end
  if ~threw
    fprintf('Expected an error but none was thrown: %s\n', label);
  end
  assert(threw, sprintf('Expected an error but none was thrown: %s', label));
end