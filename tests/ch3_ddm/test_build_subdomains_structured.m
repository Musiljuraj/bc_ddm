function test_build_subdomains_structured()
%TEST_BUILD_SUBDOMAINS_STRUCTURED  Sanity checks for structured subdomain partitioning (Chapter 3.1).

  setup_paths();
  %addpath('../../src/fem/mesh');
  %addpath('../../src/fem/bc');
  %addpath('../../src/ddm/partition');

  n = 6;          % mesh parameter (must be divisible by nSubX,nSubY)
  nSubX = 2;
  nSubY = 3;

  [p, t, bnd] = mesh_unit_square_P1(n);

  [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY);

  assert(ddm.Nsub == nSubX*nSubY);
  assert(numel(sub) == ddm.Nsub);

  % Check that each element is assigned exactly once.
  all_elems = [];
  for i = 1:ddm.Nsub
    all_elems = [all_elems; sub(i).elems(:)]; %#ok<AGROW>
  end
  all_elems_sorted = sort(all_elems);
  assert(numel(all_elems_sorted) == size(t,1));
  assert(all(all_elems_sorted(:) == (1:size(t,1)).'));

  % Check loc2glob validity and glob2loc consistency.
  for i = 1:ddm.Nsub
    g = sub(i).loc2glob(:);
    assert(all(g >= 1 & g <= ddm.nFree));
    assert(issorted(g));

    % Invertibility on the local set.
    for k = 1:numel(g)
      assert(sub(i).glob2loc(g(k)) == k);
    end
  end

  fprintf('PASS: test_build_subdomains_structured (n=%d, nSubX=%d, nSubY=%d)\n', n, nSubX, nSubY);
end
