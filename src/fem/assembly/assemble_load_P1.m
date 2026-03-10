function F = assemble_load_P1(p, t, f_handle)
%ASSEMBLE_LOAD_P1 Assemble the global P1 load vector on the full mesh.
% Thesis link: Chapter 3.3–3.4 (discrete Ritz system, assembly, boundary treatment).
% Each element contributes a local vector from `triP1_load`, which is summed
% at global node indices. The result is returned before Dirichlet elimination.
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
    error('assemble_load_P1: t must be Ne x 3.'); % FIXED (was broken across lines)
  end
  if ~isa(f_handle, 'function_handle')
    error('assemble_load_P1: f_handle must be a function handle, e.g. @(x,y) ...');
  end

  % Harden input types/values to prevent silent coercions (char/logical/NaN/Inf). % ADDED
  if ~isnumeric(p) || islogical(p) || ~isreal(p) || any(~isfinite(p(:)))         % ADDED
    error('assemble_load_P1: p must be numeric, real, finite (non-logical).');   % ADDED
  end                                                                            % ADDED
  if ~isnumeric(t) || islogical(t) || ~isreal(t) || any(~isfinite(t(:)))         % ADDED
    error('assemble_load_P1: t must be numeric, real, finite (non-logical).');   % ADDED
  end                                                                            % ADDED

  % Normalize to double for consistent downstream arithmetic/indexing.            % ADDED
  p = double(p);                                                                  % ADDED
  t = double(t);                                                                  % ADDED

  N  = size(p,1);
  Ne = size(t,1);

  % Connectivity must be valid indices 1..N and integer-valued.
  % (Use fix after finiteness checks to avoid NaN/Inf bypass.)                     % CHANGED
  if ~isempty(t)                                                                  % ADDED
    if any(t(:) < 1) || any(t(:) > N) || any(t(:) ~= fix(t(:)))                  % CHANGED
      error('assemble_load_P1: t must contain valid integer node indices in 1..N.');
    end
  end                                                                             % ADDED

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