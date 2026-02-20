function [p, t, bnd] = mesh_unit_square_P1(n)
%MESH_UNIT_SQUARE_P1  Uniform P1 triangulation of the unit square [0,1]x[0,1].
%
% Link to thesis:
%   Chapter 2.2 (Discrete approximation space on the square domain).
%   Provides the geometric mesh (nodes + triangles) on which the nodal P1 basis
%   functions (hat functions) are defined.
%
% Theoretical role:
%   The discrete space V_n is constructed on a triangulation of Ω = (0,1)^2.
%   Each global degree of freedom corresponds to one mesh node, and each
%   triangular element carries three local P1 basis functions associated with
%   its vertices.
%
%   [p, t, bnd] = mesh_unit_square_P1(n)
%
%   Input:
%   n  - number of steps (uniformly distanced nodes) per axis (integer, n >= 1)
%
% Outputs and meaning:
%   p : N×2 array of node coordinates (global DOF locations).
%   t : Ne×3 array of element (triangle) connectivities (global node indices per triangle).
%       Nodes in each triangle are ordered counter-clockwise.
%   bnd : boundary index sets used in Chapter 2.4 (Dirichlet/Neumann handling),
%         and mesh parameter h.
%          .dirichlet_nodes : column vector of global node indices with prescribed Dirichlet BC
%          .neumann_edges   : array of node index pairs (optional use) which correspond to edges with Neumann BC
%          .h               : mesh size (1/n)
%          .n               : steps per axis
%
% Implementation outline:
%   1) Create a uniform (n+1)×(n+1) grid on [0,1]×[0,1] and number nodes
%      consistently (left-to-right, bottom-to-top).
%   2) Split each grid square into two triangles with consistent orientation.
%   3) Identify boundary nodes (Dirichlet) and edges (Neumann) and return them in a small structure for
%      later restriction/elimination of Dirichlet DOFs.

  % check correct input
  if nargin ~= 1
    error('mesh_unit_square_P1: expected exactly one input argument n.');
  end
  if ~(isscalar(n) && n == floor(n) && n >= 1)
    error('mesh_unit_square_P1: n must be an integer >= 1.');
  end

  %Mesh parameters
  h = 1.0 / n;          % mesh step
  np = n + 1;           % points per axis
  N  = np * np;         % number of nodes
  Ne = 2 * n * n;       % number of triangles

  % create node coordinates p (index/number, x-coor, y-coor)
  % numbering = left-to-right, bottom-to-top
  p = zeros(N, 2);
  idx = 1;
  for j = 1:np
    y = (j-1) * h;
    for i = 1:np
      x = (i-1) * h;
      p(idx, :) = [x, y];
      idx = idx + 1;
    end
  end

  % helper function: mapping nodes i,j coordinates to node's index (used in triangle connectivity)
  id = @(i,j) i + (j-1)*np;

  % --- triangle connectivity t (counter-clockwise)
  % Each square cell (i,j) with corners:
  %   bl (bottomLeft) = (i, j), br (bottomRight) = (i+1, j), tl (topLeft) = (i, j+1), tr (topRight) = (i+1, j+1)
  % Split along diagonal bottom-right -> top-left:
  %   T1 = [bl, br, tl]
  %   T2 = [br, tr, tl]
  t = zeros(Ne, 3);
  e = 1;
  for j = 1:n
    for i = 1:n
      bl = id(i,   j);
      br = id(i+1, j);
      tl = id(i,   j+1);
      tr = id(i+1, j+1);

      t(e, :) = [bl, br, tl]; e = e + 1;  % counter-clockwise
      t(e, :) = [br, tr, tl]; e = e + 1;  % counter-clockwise
    end
  end

  % --- boundary marking
  % Thesis model is square domain with Dirichlet boundaries on left and right edges (x=0 and x=1)
  % and Neumann boundarier on bottom and top edges (y=0 and y=1) - stored as edges list
  tol = 1e-12; 

  x = p(:,1); y = p(:,2);

  %testing whether the node is actually on boundary using tollerance (due to possible rounding errors)
  dirichlet_nodes = find( (abs(x - 0.0) < tol) | (abs(x - 1.0) < tol) ); 
  dirichlet_nodes = unique(dirichlet_nodes(:)); %actually not needed for my specific thesis setup, but better be safe than sorry :-)


  % Neumann edges (bottom and top boundary) as pairs of node indices 
  % (Neumann boundary conditions are stored as edges, not nodes)
  % For my specific thesis case (homogenous Neumann BC actually not needed, but implemented anyway)
  % bottom: j=1, i=1..n -> edges (i,1)-(i+1,1)
  % top:    j=np, i=1..n -> edges (i,np)-(i+1,np)
  neumann_edges = zeros(2*n, 2);
  k = 1;
  % bottom edges
  j = 1;
  for i = 1:n
    neumann_edges(k,:) = [id(i,j), id(i+1,j)];
    k = k + 1;
  end
  % top edges
  j = np;
  for i = 1:n
    neumann_edges(k,:) = [id(i,j), id(i+1,j)];
    k = k + 1;
  end

  %packing all relevant boundary related data
  bnd = struct();
  bnd.dirichlet_nodes = dirichlet_nodes;
  bnd.neumann_edges   = neumann_edges;
  bnd.h               = h;
  bnd.n               = n;
end
