function K = assemble_stiffness_P1(p, t)
%ASSEMBLE_STIFFNESS_P1  Assemble the global P1 stiffness matrix on a triangulation.
%
% Link to thesis:
%   Chapter 2.3.3 (Stiffness matrix assembly) and the element-wise representation
%   K_ij = Σ_e ∫_{T_e} ∇φ_i·∇φ_j dx.
%
% Theoretical role:
%   This function constructs the global stiffness matrix K that appears in the
%   discrete Ritz / FEM system
%   The global stiffness matrix K represents the bilinear form A(·,·) restricted
%   to the discrete space V_n. Each element contributes a 3×3 matrix K^(e), and
%   shared nodes imply that global entries receive summed contributions from
%   multiple elements.
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

% IMPORTANT NOTE: SPARSE ASSEMBLY VIA TRIPLETS (I,J,V)
% --------------------------------------------------------------
% The stiffness matrix K is sparse because each basis function phi_i has local
% support: it only overlaps with basis functions associated with neighboring
% nodes. Therefore, each row of K has only O(1) nonzeros independent of N.
%
% Building sparse matrices efficiently in Octave/Matlab is done by providing
% triplets (I,J,V) and calling:
%
%       K = sparse(I, J, V, N, N)
%
% Semantics:
%   For each k, the value V(k) is ADDED to K(I(k), J(k)).
%   If indices repeat (which they will in FEM assembly), sparse() sums them.
%
% -------------------------------------------------------------------------

  % -----------------------------
  % 1) Basic input validation
  % -----------------------------
  if nargin ~= 2
    error('assemble_stiffness_P1: expected inputs (p, t).');
  end
  if ~ismatrix(p) || size(p,2) ~= 2
    error('assemble_stiffness_P1: p must be N x 2.');
  end
  if ~ismatrix(t) || size(t,2) ~= 3
    error('assemble_stiffness_P1: t must be Ne x 3.');
  end

  N  = size(p,1);    % total number of global degrees of freedom (nodes)
  Ne = size(t,1);    % number of triangular elements

  % Connectivity must consist of valid global node indices in 1..N.
  if any(t(:) < 1) || any(t(:) > N) || any(abs(t(:) - round(t(:))) > 0)
    error('assemble_stiffness_P1: t must contain valid integer node indices in 1..N.');
  end

  % -------------------------------------------------------------
  % 2) Allocate triplet storage (I,J,V) for sparse matrix assembly
  % -------------------------------------------------------------
  % Each triangle contributes a 3×3 local matrix Ke -> 9 scalar entries.
  %
  % We will store these 9 entries as 9 triplets:
  %    (global_row, global_col, value)
  %
  % Over all Ne elements we therefore store exactly 9*Ne triplets.
  %
  % After the loop, sparse(I,J,V,N,N) will automatically sum repeated (I,J).
  %
  I = zeros(9*Ne, 1);   % global row indices of contributions
  J = zeros(9*Ne, 1);   % global column indices of contributions
  V = zeros(9*Ne, 1);   % contribution values (Ke entries)

  % idx points to the next free block of 9 positions in (I,J,V).
  idx = 1;

  % -----------------------------
  % 3) Element loop (assembly)
  % -----------------------------
  for e = 1:Ne

    % -----------------------------------------------------------
    % 3.1) Local-to-global mapping for this triangle
    % -----------------------------------------------------------
    % nodes are the GLOBAL node numbers of the triangle vertices:
    %   nodes = [i1, i2, i3]
    %
    % In FEM theory, this represents the mapping:
    %   local DoF 1 -> global DoF i1
    %   local DoF 2 -> global DoF i2
    %   local DoF 3 -> global DoF i3
    %
    nodes = t(e,:);          % 1×3 vector of global indices

    % xy are the coordinates of the triangle vertices in the same order
    % as 'nodes'. This is the geometric data needed by the element routine.
    xy = p(nodes,:);         % 3×2

    % -----------------------------------------------------------
    % 3.2) Compute local stiffness matrix Ke for this triangle
    % -----------------------------------------------------------
    % Ke(m,n) = ∫_{T_e} ∇phi_m · ∇phi_n dx,  m,n = 1..3 (local basis functions)
    %
    % triP1_stiffness(xy) evaluates this exactly for P1 triangles.
    %
    Ke = triP1_stiffness(xy);    % 3×3

    % -----------------------------------------------------------
    % 3.3) Convert the 3×3 local matrix into 9 global triplets
    % -----------------------------------------------------------
    % Theory: 
    %   Add Ke(m,n) into global K at position:
    %       K( nodes(m), nodes(n) ) += Ke(m,n)
    %
    % That means we need all pairs (nodes(m), nodes(n)) for m,n = 1..3.
    %
    % We generate these 9 index pairs using meshgrid(nodes, nodes) which produces two 3×3 matrices ii and jj such that:
    %
    %   ii =
    %     i1  i2  i3
    %     i1  i2  i3
    %     i1  i2  i3
    %
    %   jj =
    %     i1  i1  i1
    %     i2  i2  i2
    %     i3  i3  i3
    %
    % The pairs (ii(m,n), jj(m,n)) then enumerate all (nodes(column), nodes(row))
    % combinations. Then ii(:), jj(:) flatten these 3×3 matrices into 9×1
    % vectors (column-wise).
    %
    % IMPORTANT:
    % We must use the SAME flattening order for indices and Ke(:) values.
    % Since Ke(:) also flattens column-wise, the entries match correctly:
    %   Ke(:) = [Ke(1,1); Ke(2,1); Ke(3,1); Ke(1,2); ... ; Ke(3,3)]
    %
    % The corresponding indices (I,J) must follow the same order.
    %
    [ii, jj] = meshgrid(nodes, nodes);

    % Store the 9 contributions of this element into the next block of (I,J,V).
    I(idx:idx+8) = ii(:);     % 9×1 global row indices
    J(idx:idx+8) = jj(:);     % 9×1 global col indices
    V(idx:idx+8) = Ke(:);     % 9×1 values

    % Advance pointer by 9 for the next element.
    idx = idx + 9;
  end

  % -------------------------------------------------------------
  % 4) Build sparse global stiffness matrix K from triplets
  % -------------------------------------------------------------
  % This call performs:
  %   K(I(k), J(k)) += V(k)   for k = 1..9*Ne
  %
  % If multiple elements contribute to the same global entry (i,j),
  % sparse() sums them automatically.
  %
  K = sparse(I, J, V, N, N);

end
