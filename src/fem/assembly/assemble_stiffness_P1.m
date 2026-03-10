function K = assemble_stiffness_P1(p, t)
%ASSEMBLE_STIFFNESS_P1 Assemble the global P1 stiffness matrix on a triangulation.
% Thesis link: Chapter 3.3 (element stiffness and global assembly).
% The matrix is built from local `triP1_stiffness` contributions using
% sparse triplet assembly and is returned before Dirichlet elimination.
%
% Inputs:
%   p : N×2 of global node coordinates.
%   t : Ne×3 array of triangle connectivities (global node indices).
%     (the ordering is CCW in my mesh generator, but Ke is invariant)
% Outputs:
%   K : N×N, final sparse global stiffness matrix (but before removing Dirichlet nodes)
%
% Implementation outline:
%   1) For each element e, extract its vertex coordinates from p and t.
%   2) Compute K^(e) using triP1_stiffness.
%   3) Map local indices {1,2,3} to global node indices t(e,:) and add all
%      9 contributions into the global sparse matrix.
%   The assembly realizes the local-to-global summation described in the thesis.

  % -----------------------------
  % 1) Basic input validation
  % -----------------------------
  if nargin ~= 2
    error('assemble_stiffness_P1: expected inputs (p, t).');
  end

  % ADDED: stronger type/real/finite checks (avoid silent ASCII/NaN/Inf issues)
  if ~isnumeric(p) || ~isreal(p)
    error('assemble_stiffness_P1: p must be a real numeric array.');
  end
  if ~ismatrix(p) || size(p,2) ~= 2
    error('assemble_stiffness_P1: p must be N x 2.');
  end
  if any(~isfinite(p(:)))
    error('assemble_stiffness_P1: p must contain only finite values.');
  end

  % ADDED: stronger checks for t (numeric/real/finite) before index validation
  if ~isnumeric(t) || ~isreal(t)
    error('assemble_stiffness_P1: t must be a real numeric array.');
  end
  if ~ismatrix(t) || size(t,2) ~= 3
    error('assemble_stiffness_P1: t must be Ne x 3.');
  end
  if any(~isfinite(t(:)))
    error('assemble_stiffness_P1: t must contain only finite values.');
  end

  N  = size(p,1);    % total number of global degrees of freedom (nodes)
  Ne = size(t,1);    % number of triangular elements

  % Connectivity must consist of valid global node indices in 1..N.
  % CHANGED: keep the original intent but use a clearer integer check
  if any(t(:) ~= round(t(:))) || any(t(:) < 1) || any(t(:) > N)
    error('assemble_stiffness_P1: t must contain valid integer node indices in 1..N.');
  end

  % ADDED: reject degenerate connectivity with repeated vertices (common silent failure mode)
  if Ne > 0
    if any(t(:,1) == t(:,2) | t(:,1) == t(:,3) | t(:,2) == t(:,3))
      error('assemble_stiffness_P1: each triangle must reference 3 distinct node indices.');
    end
  end

  % -------------------------------------------------------------
  % 2) Allocate triplet storage (I,J,V) for sparse matrix assembly
  % -------------------------------------------------------------
  % Each triangle contributes a 3×3 local matrix Ke -> 9 scalar entries.
  I = zeros(9*Ne, 1);
  J = zeros(9*Ne, 1);
  V = zeros(9*Ne, 1);

  idx = 1;

  % -----------------------------
  % 3) Element loop (assembly)
  % -----------------------------
  for e = 1:Ne

    nodes = t(e,:);          % 1×3 vector of global indices
    xy = p(nodes,:);         % 3×2

    Ke = triP1_stiffness(xy);    % 3×3

    % ADDED: sanity check the element routine output (catches accidental API drift)
    if ~ismatrix(Ke) || any(size(Ke) ~= [3 3]) || ~isnumeric(Ke) || ~isreal(Ke)
      error('assemble_stiffness_P1: triP1_stiffness must return a real numeric 3x3 matrix.');
    end
    if any(~isfinite(Ke(:)))
      error('assemble_stiffness_P1: triP1_stiffness returned NaN/Inf (degenerate element?).');
    end

    % -----------------------------------------------------------
    % 3.3) Convert the 3×3 local matrix into 9 global triplets
    % -----------------------------------------------------------
    % Theory:
    %   K(nodes(m), nodes(n)) += Ke(m,n),  m,n = 1..3
    %
    % FIXED: use ndgrid (not meshgrid) so that (I,J) flattening matches Ke(:)
    %        in column-major order:
    %          Ke(:) = [Ke(1,1); Ke(2,1); Ke(3,1); Ke(1,2); ... ; Ke(3,3)]
    %        With ndgrid(nodes,nodes):
    %          ii(:) = [nodes(1); nodes(2); nodes(3); nodes(1); ...]  (row indices)
    %          jj(:) = [nodes(1); nodes(1); nodes(1); nodes(2); ...]  (col indices)
    %
    %        This prevents accidentally assembling Ke' (transpose) due to an
    %        index/value ordering mismatch.
    [ii, jj] = ndgrid(nodes, nodes);  % FIXED

    I(idx:idx+8) = ii(:);
    J(idx:idx+8) = jj(:);
    V(idx:idx+8) = Ke(:);

    idx = idx + 9;
  end

  % -------------------------------------------------------------
  % 4) Build sparse global stiffness matrix K from triplets
  % -------------------------------------------------------------
  K = sparse(I, J, V, N, N);

end