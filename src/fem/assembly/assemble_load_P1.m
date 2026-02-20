function F = assemble_load_P1(p, t, f_handle)
%ASSEMBLE_LOAD_P1  Assemble global load vector for P1 FEM on triangles.
%
% Link to thesis:
%   Chapter 2.3.1 (load vector entries) together with the element-wise sum
%   f_i = Σ_e ∫_{T_e} φ_i f dx and the assembly concept from Chapter 2.3.3.
%
% Theoretical role:
%   The global vector F represents the linear functional ℓ(v)=∫_Ω v f dx evaluated
%   at basis functions (for homogeneous Neumann BCs: no boundary term). 
%   Each triangle contributes a local vector f^(e) that is
%   added to F at the global indices of the triangle vertices.
%
% Input:
%   p        : N×2 node coordinates (global).
%   t        : Ne×3 triangle connectivities (global node indices).
%   f_handle : function handle for the load (force) term f(x,y).
% Outputs:
%   F        : N×1 assembled load vector on all nodes (prior to Dirichlet BC reduction).
%
% Implementation outline:
%   1) For each element e, extract vertex coordinates.
%   2) Compute the local vector f^(e) using triP1_load (centroid quadrature).
%   3) Add the three local entries into F at indices t(e,:).

% IMPORTANT NOTE: SPARSE ASSEMBLY VIA sparse constructon
% --------------------------------------------------------------
% For stiffness matrix assembly we  used sparse triplets (I,J,V).
% For the load vector (a column vector), we can use an analogous sparse
% constructor:
%
%   F = sparse(I, 1, V, N, 1)
%
% where:
%   - I stores global row indices (node numbers)
%   - column index is always 1 (because F is N×1)
%   - V stores the values to be summed at those rows
%
% Repeated indices in I are automatically summed, exactly as required by FEM.

  % -----------------------------
  % 1) Input validation
  % -----------------------------
  if nargin ~= 3
    error('assemble_load_P1: expected inputs (p, t, f_handle).');
  end
  if ~ismatrix(p) || size(p,2) ~= 2
    error('assemble_load_P1: p must be N x 2.');
  end
  if ~ismatrix(t) || size(t,2) ~= 3
    error('assemble_load_P1: t must be Ne x 3.');
  end
  if ~isa(f_handle, 'function_handle')
    error('assemble_load_P1: f_handle must be a function handle, e.g. @(x,y) ...');
  end

  N  = size(p,1);
  Ne = size(t,1);

  % Connectivity must be valid indices 1..N and integer-valued.
  if any(t(:) < 1) || any(t(:) > N) || any(abs(t(:) - round(t(:))) > 0)
    error('assemble_load_P1: t must contain valid integer node indices in 1..N.');
  end

  % ---------------------------------------------------------
  % 2) Preallocate triplets for sparse vector assembly
  % ---------------------------------------------------------
  % Each element contributes 3 entries (one per local basis function).
  %
  % I(k): row index (global node index)
  % V(k): value to add to F(I(k))
  %
  I = zeros(3*Ne, 1);
  V = zeros(3*Ne, 1);

  idx = 1;  % pointer to next free block of length 3

  % -----------------------------
  % 3) Element loop (assembly)
  % -----------------------------
  for e = 1:Ne

    % Global node numbers of element e (local-to-global mapping)
    nodes = t(e,:);      % 1×3

    % Vertex coordinates of the triangle
    xy = p(nodes,:);     % 3×2

    % Local element load vector fe (3×1)
    fe = triP1_load(xy, f_handle);

    % Store triplets for this element:
    % - rows are the global node indices
    % - values are fe entries
    I(idx:idx+2) = nodes(:);   % nodes(:) converts 1×3 into 3×1
    V(idx:idx+2) = fe(:);      % fe is already 3×1, but fe(:) is explicit

    idx = idx + 3;
  end

  % ---------------------------------------------------------
  % 4) Build global load vector using sparse summation
  % ---------------------------------------------------------
  % sparse(I, 1, V, N, 1) creates an N×1 sparse vector where:
  %   F(I(k), 1) += V(k).
  %
  % This automatically sums contributions for nodes shared by multiple
  % triangles, exactly matching the FEM formula F_i = Σ_e ∫_{T_e} φ_i f.
  %
  F = sparse(I, ones(size(I)), V, N, 1);

end
