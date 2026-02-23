function test_fem_block_integration()
%TEST_FEM_BLOCK_INTEGRATION  Intermediate/global test of the FEM block.
%
% This test exercises the full FEM chain:
%   mesh -> element routines -> assembly -> Dirichlet elimination -> solve
%
% It includes:
%   (1) Manufactured-solution refinement ladder (n = 8,16,32) to verify
%       meaningful PDE results and convergence.
%   (2) Structural invariants (symmetry, SPD of reduced matrix, nullspace of K).
%   (3) Negative/fault-injection cases to ensure robust failure modes.
%
% Octave-compatible: plain assertions, no matlab.unittest.

  % Keep tests runnable standalone (harmless if caller already set paths).
  if exist('setup_paths','file') == 2
    setup_paths();
  end

  % Manufactured solution compatible with BC choice in mesh_unit_square_P1:
  % Dirichlet on x=0 and x=1, homogeneous Neumann on y=0 and y=1.
  u_exact = @(x,y) sin(pi*x) .* cos(pi*y);
  f_rhs   = @(x,y) 2*pi^2 * sin(pi*x) .* cos(pi*y);

  n_list = [8 16 32];
  outs = cell(numel(n_list), 1);
  err = zeros(numel(n_list), 1);

  % --- Positive path: solve on refinement ladder ---
  for k = 1:numel(n_list)
    cfg = struct();
    cfg.n = n_list(k);
    cfg.f_handle = f_rhs;
    cfg.u_exact_handle = u_exact;
    cfg.solve = true;

    out = run_fem_block_case(cfg);
    outs{k} = out;

    d = out.diagnostics;
    err(k) = d.errL2_rel;

    % Basic invariants
    assert(d.rel_sym_K   < 1e-12, sprintf('K not symmetric enough for n=%d (rel=%g).', cfg.n, d.rel_sym_K));
    assert(d.rel_sym_Kff < 1e-12, sprintf('Kff not symmetric enough for n=%d (rel=%g).', cfg.n, d.rel_sym_Kff));
    assert(d.chol_ok, sprintf('Reduced matrix Kff not SPD (chol failed) for n=%d.', cfg.n));

    % K*ones ~= 0 (nullspace) check (scale-aware, but not too strict)
    assert(d.rel_null_K < 1e-10, sprintf('K*ones nullspace check failed for n=%d (rel=%g).', cfg.n, d.rel_null_K));

    % Dirichlet nodes should be exactly zero in reconstructed solution
    dn = out.bnd.dirichlet_nodes(:);
    assert(norm(out.u_full(dn), inf) == 0, sprintf('Dirichlet nodes not exactly zero for n=%d.', cfg.n));

    % Residual should be small for direct solve
    assert(d.rel_resid < 1e-12, sprintf('Linear solve residual too large for n=%d (rel=%g).', cfg.n, d.rel_resid));

    % Simple physical sanity: bottom and top center should have opposite signs
    % for u=sin(pi x)cos(pi y): at y=0 positive, at y=1 negative.
    if isfield(d,'probe') && isfield(d.probe,'idx_bot') && ~isnan(d.probe.idx_bot)
      assert(d.probe.u_bot > 0, sprintf('Unexpected sign at bottom center for n=%d (u=%g).', cfg.n, d.probe.u_bot));
      assert(d.probe.u_top < 0, sprintf('Unexpected sign at top center for n=%d (u=%g).', cfg.n, d.probe.u_top));
    end
  end

  % Convergence checks (conservative): error should decrease with refinement.
  assert(err(2) < 0.95*err(1), sprintf('Expected error decrease from n=%d to n=%d (err %g -> %g).', n_list(1), n_list(2), err(1), err(2)));
  assert(err(3) < 0.95*err(2), sprintf('Expected error decrease from n=%d to n=%d (err %g -> %g).', n_list(2), n_list(3), err(2), err(3)));

  % Self-convergence check: u_n should approach u_2n when restricting fine solution
  % onto the coarse grid nodes (nested structured grid).
  diff12 = rel_diff_coarse_vs_fine(outs{1}, outs{2});
  diff23 = rel_diff_coarse_vs_fine(outs{2}, outs{3});
  assert(diff23 < 0.90*diff12, sprintf('Self-convergence not improving enough (diff12=%g, diff23=%g).', diff12, diff23));

  % --- Negative path / fault injection ---
  negative_cases(f_rhs);

  fprintf('PASS: test_fem_block_integration\n');
end

% ========================================================================
% Helpers
% ========================================================================

function rel = rel_diff_coarse_vs_fine(out_c, out_f)
  % Compare u on coarse mesh vs restricted u on fine mesh for free nodes.
  nc = out_c.bnd.n;
  nf = out_f.bnd.n;
  assert(nf == 2*nc, 'Expected nested refinement nf = 2*nc.');

  map = coarse_to_fine_map(nc);

  free_c = out_c.free_nodes(:);
  uc = out_c.u_full(free_c);
  uf = out_f.u_full(map(free_c));

  rel = norm(uc - uf, 2) / max(1e-16, norm(uf, 2));
end

function map = coarse_to_fine_map(nc)
  % Map coarse node index -> fine node index for structured nested meshes.
  % Coarse has np_c = nc+1 nodes/axis, fine has np_f = 2*nc+1.
  np_c = nc + 1;
  np_f = 2*nc + 1;
  N_c  = np_c * np_c;

  map = zeros(N_c, 1);
  for k = 1:N_c
    i = mod(k-1, np_c) + 1;        % 1..np_c
    j = floor((k-1) / np_c) + 1;  % 1..np_c

    i_f = 2*(i-1) + 1;
    j_f = 2*(j-1) + 1;

    map(k) = i_f + (j_f-1)*np_f;
  end
end

function negative_cases(f_rhs)
  % A small set of end-to-end failure-mode checks.

  % 1) Degenerate triangle should trigger triP1_stiffness:degenerate during assembly.
  cfg = struct();
  cfg.n = 8;
  cfg.f_handle = f_rhs;
  cfg.solve = false;
  cfg.mutate = @mutate_make_degenerate_triangle;
  assert_throws(@() run_fem_block_case(cfg), {'triP1_stiffness:degenerate', 'degenerate triangle'});

  % 2) Out-of-range connectivity should be rejected by assemble_stiffness_P1.
  cfg = struct();
  cfg.n = 8;
  cfg.f_handle = f_rhs;
  cfg.solve = false;
  cfg.mutate = @mutate_bad_connectivity_zero;
  assert_throws(@() run_fem_block_case(cfg), {'valid integer node indices', 't must contain valid'});

  % 3) Repeated vertex in a triangle should be rejected.
  cfg = struct();
  cfg.n = 8;
  cfg.f_handle = f_rhs;
  cfg.solve = false;
  cfg.mutate = @mutate_repeated_vertex;
  assert_throws(@() run_fem_block_case(cfg), {'distinct node indices', '3 distinct'});

  % 4) NaN RHS should be rejected by triP1_load.
  cfg = struct();
  cfg.n = 8;
  cfg.f_handle = @(x,y) NaN;
  cfg.solve = false;
  assert_throws(@() run_fem_block_case(cfg), {'finite scalar', 'f_handle must return'});

  % 5) Elimination with all nodes Dirichlet should fail.
  cfg = struct();
  cfg.n = 8;
  cfg.f_handle = f_rhs;
  cfg.solve = false;
  cfg.mutate = @mutate_all_dirichlet;
  assert_throws(@() run_fem_block_case(cfg), {'no free degrees of freedom'});
end

function [p, t, bnd] = mutate_make_degenerate_triangle(p, t, bnd)
  % Force the first triangle to be degenerate by moving one vertex onto the
  % line segment between the other two.
  tri = t(1,:);
  a = tri(1);
  b = tri(2);
  c = tri(3);
  p(c,:) = 0.5*(p(a,:) + p(b,:));
end

function [p, t, bnd] = mutate_bad_connectivity_zero(p, t, bnd)
  t(1,1) = 0;
end

function [p, t, bnd] = mutate_repeated_vertex(p, t, bnd)
  t(1,2) = t(1,1);
end

function [p, t, bnd] = mutate_all_dirichlet(p, t, bnd)
  bnd.dirichlet_nodes = (1:size(p,1)).';
end

function assert_throws(fun, expected)
  % Assert that fun throws, and that identifier/message matches expected.
  did_throw = false;
  try
    fun();
  catch err
    did_throw = true;

    if nargin >= 2 && ~isempty(expected)
      if ischar(expected)
        expected = {expected};
      end
      ok = false;
      for k = 1:numel(expected)
        pat = expected{k};
        if isfield(err,'identifier') && ~isempty(err.identifier) && ~isempty(strfind(err.identifier, pat))
          ok = true; break;
        end
        if isfield(err,'message') && ~isempty(strfind(err.message, pat))
          ok = true; break;
        end
      end
      assert(ok, sprintf('Error did not match expected patterns. Got id="%s", msg="%s"', safe_str(err,'identifier'), safe_str(err,'message')));
    end
  end
  assert(did_throw, 'Expected an error, but no error was thrown.');
end

function s = safe_str(err, field)
  s = '';
  if isfield(err, field)
    try
      s = err.(field);
    catch
      s = '';
    end
  end
  if isempty(s)
    s = '';
  end
end