function test_mesh_unit_square_P1()
%TEST_MESH_UNIT_SQUARE_P1  Unit tests for uniform P1 unit-square mesh.

  setup_paths();

  % ---- Signature consistency
  assert(nargin('mesh_unit_square_P1') == 1, 'mesh_unit_square_P1 should take exactly 1 input.');
  assert(nargout('mesh_unit_square_P1') == 3, 'mesh_unit_square_P1 should define exactly 3 outputs.');

  % ---- Input validation
  % For wrong number of arguments, Octave may throw its own message.
  assert_throws(@() mesh_unit_square_P1(), '');
  assert_throws(@() mesh_unit_square_P1(1, 2), '');

  % Your function's own message contains "integer >= 1" so we can check that substring.
  assert_throws(@() mesh_unit_square_P1(0),     'integer >= 1');
  assert_throws(@() mesh_unit_square_P1(-3),    'integer >= 1');
  assert_throws(@() mesh_unit_square_P1(2.5),   'integer >= 1');
  assert_throws(@() mesh_unit_square_P1([2 3]), 'integer >= 1');

  % These should now pass because your mesh_unit_square_P1.m rejects char/logical.
  assert_throws(@() mesh_unit_square_P1('4'),   'integer >= 1');
  assert_throws(@() mesh_unit_square_P1(true),  'integer >= 1');

  % Inf should error (message may vary across Octave versions, so don't enforce substring).
  assert_throws(@() mesh_unit_square_P1(Inf), '');

  % ---- Functional tests for several n
  ns = [1, 2, 4];
  for kk = 1:numel(ns)
    n = ns(kk);
    [p, t, bnd] = mesh_unit_square_P1(n);

    np = n + 1;
    N  = np*np;
    Ne = 2*n*n;
    h  = 1.0/n;

    % Basic shapes
    assert(isequal(size(p), [N, 2]),  sprintf('p has wrong size for n=%d', n));
    assert(isequal(size(t), [Ne, 3]), sprintf('t has wrong size for n=%d', n));
    assert(isstruct(bnd),             sprintf('bnd should be a struct for n=%d', n));

    % bnd fields and values
    assert(isfield(bnd,'dirichlet_nodes'), sprintf('bnd.dirichlet_nodes missing for n=%d', n));
    assert(isfield(bnd,'neumann_edges'),   sprintf('bnd.neumann_edges missing for n=%d', n));
    assert(isfield(bnd,'h'),               sprintf('bnd.h missing for n=%d', n));
    assert(isfield(bnd,'n'),               sprintf('bnd.n missing for n=%d', n));
    assert_close(bnd.h, h, 1e-14, 'bnd.h mismatch');
    assert(bnd.n == n, sprintf('bnd.n mismatch for n=%d', n));

    % ---- Exact grid coordinates and ordering
    tol = 5e-13;
    p_expected = zeros(N, 2);
    idx = 1;
    for j = 0:n
      y = j*h;
      for i = 0:n
        x = i*h;
        p_expected(idx,:) = [x, y];
        idx = idx + 1;
      end
    end
    assert(max(abs(p(:) - p_expected(:))) < tol, sprintf('p differs from expected grid/order for n=%d', n));

    % Corners
    assert_close(p(1,:),      [0,0], tol, 'bottom-left coordinate wrong');
    assert_close(p(np,:),     [1,0], tol, 'bottom-right coordinate wrong');
    assert_close(p(N-np+1,:), [0,1], tol, 'top-left coordinate wrong');
    assert_close(p(N,:),      [1,1], tol, 'top-right coordinate wrong');

    % ---- Triangle connectivity sanity
    assert(all(t(:) >= 1) && all(t(:) <= N), sprintf('Triangle connectivity out of range for n=%d', n));

    % Each triangle must have 3 distinct vertices
    assert(all(t(:,1) ~= t(:,2) & t(:,1) ~= t(:,3) & t(:,2) ~= t(:,3)), ...
           sprintf('Triangle with repeated vertex detected for n=%d', n));

    % Orientation and area checks
    p1 = p(t(:,1),:); p2 = p(t(:,2),:); p3 = p(t(:,3),:);
    signed_area = 0.5 * ( (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) - (p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1)) );

    assert(all(signed_area > 0), sprintf('Non-CCW (or degenerate) triangle detected for n=%d', n));

    area_expected = (h*h)/2;
    assert(max(abs(signed_area - area_expected)) < 5e-13, sprintf('Triangle area mismatch for n=%d', n));

    total_area = sum(signed_area);
    assert_close(total_area, 1.0, 5e-12, 'Total mesh area should be 1');

    % ---- Exact expected triangle set (order-independent)
    tri_can = sort(t, 2);
    tri_can = sortrows(tri_can);

    tri_exp = zeros(Ne, 3);
    e = 1;
    for j = 1:n
      for i = 1:n
        bl = id_ij(i,   j,   np);
        br = id_ij(i+1, j,   np);
        tl = id_ij(i,   j+1, np);
        tr = id_ij(i+1, j+1, np);
        tri_exp(e,:) = sort([bl, br, tl]); e = e + 1;
        tri_exp(e,:) = sort([br, tr, tl]); e = e + 1;
      end
    end
    tri_exp = sortrows(tri_exp);

    assert(isequal(tri_can, tri_exp), sprintf('Triangle set mismatch for n=%d', n));

    % ---- Boundary: Dirichlet nodes (left/right), exact indices
    left  = (0:np-1)'*np + 1;
    right = (1:np)'*np;
    dir_exp = sort([left; right]);

    dir_act = sort(unique(bnd.dirichlet_nodes(:)));
    assert(isequal(dir_act, dir_exp), sprintf('Dirichlet node set mismatch for n=%d', n));

    % ---- Boundary: Neumann edges (bottom/top), order- and orientation-independent
    ne_act = bnd.neumann_edges;
    assert(isequal(size(ne_act), [2*n, 2]), sprintf('Neumann edges size mismatch for n=%d', n));

    bottom = [(1:n)' (2:n+1)'];
    top0 = (np-1)*np + 1;
    top = [ (top0:(top0+n-1))' (top0+1:(top0+n))' ];
    ne_exp = [bottom; top];

    assert(isequal(canon_edges(ne_act), canon_edges(ne_exp)), sprintf('Neumann edges set mismatch for n=%d', n));
  end

  fprintf('test_mesh_unit_square_P1: OK (tested n = [%s])\n', num2str(ns));
end

% ---- helpers

function idx = id_ij(i, j, np)
  idx = i + (j-1)*np;
end

function E = canon_edges(E)
  E = sort(E, 2);
  E = sortrows(E);
end

function assert_close(a, b, tol, msg)
  if nargin < 4, msg = 'Values not close'; end
  assert(max(abs(a(:) - b(:))) <= tol, msg);
end

function assert_throws(fhandle, must_contain)
  did_throw = false;
  try
    fhandle();
  catch err
    did_throw = true;
    if nargin >= 2 && ~isempty(must_contain)
      assert(~isempty(strfind(err.message, must_contain)), ...
             'Exception message did not contain expected substring.');
    end
  end
  assert(did_throw, 'Expected an exception, but none was thrown.');
end