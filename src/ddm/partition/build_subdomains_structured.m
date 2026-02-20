function [sub, ddm] = build_subdomains_structured(p, t, bnd, nSubX, nSubY)
%BUILD_SUBDOMAINS_STRUCTURED  Partition a structured unit-square mesh into rectangular subdomains.
%
% Link to thesis:
%   Chapter 3.1.1–3.1.2 (domain decomposition and induced local spaces),
%   equations (3.1), (3.6)–(3.8).
%
% This routine performs only geometric/topological bookkeeping:
%   - Partition the global element set into subdomains (non-overlapping by elements).
%   - Construct, for each subdomain, its element index set and its local DOF set.
%   - Construct local-to-global mapping loc2glob for *free* (non-Dirichlet) DOFs.
%
% Important:
%   Dirichlet nodes (as given by bnd.dirichlet_nodes) are eliminated from the
%   local DOF sets. The resulting DOFs correspond exactly to the "free" unknowns
%   in the reduced FEM system of Chapter 2.4.
%
% Inputs:
%   p     : Np x 2 node coordinates.
%   t     : Nt x 3 triangle connectivity (global node indices).
%   bnd   : boundary struct from mesh generator (must contain .dirichlet_nodes and .n).
%   nSubX : number of subdomains along x-direction.
%   nSubY : number of subdomains along y-direction.
%
% Outputs:
%   sub : struct array of length Nsub = nSubX*nSubY with fields:
%       .id        : subdomain id in 1..Nsub
%       .ix,.iy    : subdomain grid coordinates (1..nSubX, 1..nSubY)
%       .bbox      : [x0 x1 y0 y1]
%       .elems     : triangle indices belonging to the subdomain
%       .nodes     : unique global node indices used by these triangles
%       .loc2glob  : global *free-DOF* indices present in the subdomain (sorted)
%       .glob2loc  : map from global free DOF -> local index (0 if absent)
%       .nloc      : number of local free DOFs
%
%   ddm : global bookkeeping struct:
%       .nSubX, .nSubY, .Nsub
%       .nNodes, .nElems
%       .dirichlet_nodes, .free_nodes
%       .node2dof : map global node -> global free DOF id (0 for Dirichlet)
%       .dof2node : inverse map (global free DOF -> global node index)
%       .nFree    : number of free DOFs

  % -----------------------------
  % 1) Input validation
  % -----------------------------
  if nargin ~= 5
    error('build_subdomains_structured: expected inputs (p,t,bnd,nSubX,nSubY).');
  end
  if ~ismatrix(p) || size(p,2) ~= 2
    error('build_subdomains_structured: p must be Np x 2.');
  end
  if ~ismatrix(t) || size(t,2) ~= 3
    error('build_subdomains_structured: t must be Nt x 3.');
  end
  if ~isstruct(bnd) || ~isfield(bnd,'dirichlet_nodes') || ~isfield(bnd,'n')
    error('build_subdomains_structured: bnd must contain fields dirichlet_nodes and n.');
  end
  if ~(isscalar(nSubX) && nSubX == round(nSubX) && nSubX >= 1)
    error('build_subdomains_structured: nSubX must be a positive integer.');
  end
  if ~(isscalar(nSubY) && nSubY == round(nSubY) && nSubY >= 1)
    error('build_subdomains_structured: nSubY must be a positive integer.');
  end

  n = bnd.n;
  if mod(n, nSubX) ~= 0 || mod(n, nSubY) ~= 0
    error('build_subdomains_structured: mesh parameter n must be divisible by (nSubX,nSubY) for an aligned structured partition.');
  end

  nNodes = size(p,1);
  nElems = size(t,1);

  % -----------------------------
  % 2) Global DOF numbering after Dirichlet elimination (Chapter 2.4)
  %  mapping from "Dirichlet-included" geomotry to "free(non-Dirichlet)" geometry
  % -----------------------------
  dir_nodes = unique(bnd.dirichlet_nodes(:));
  all_nodes = (1:nNodes).';
  free_nodes = setdiff(all_nodes, dir_nodes);

  node2dof = zeros(nNodes,1);
  node2dof(free_nodes) = (1:numel(free_nodes)).';
  dof2node = free_nodes;

  %packing all subdomain-decomposition relevant informations into one structure
  ddm = struct();
  ddm.nSubX = nSubX;
  ddm.nSubY = nSubY;
  ddm.Nsub  = nSubX * nSubY;
  ddm.nNodes = nNodes;
  ddm.nElems = nElems;
  ddm.dirichlet_nodes = dir_nodes;
  ddm.free_nodes = free_nodes;
  ddm.node2dof = node2dof;
  ddm.dof2node = dof2node;
  ddm.nFree = numel(free_nodes);

  % -----------------------------
  % 3) Assign each element to a subdomain (by centroid)
  % -----------------------------
  % Subdomain boxes are aligned with [0,1]x[0,1]:
  %   Ω_{ix,iy} = ((ix-1)/nSubX, ix/nSubX) x ((iy-1)/nSubY, iy/nSubY)
  %
  % Elements are assigned by the location of their geometric centroid.
  %
  tol = 1e-12;

  % Precompute triangle centroids.
  nodes1 = t(:,1); nodes2 = t(:,2); nodes3 = t(:,3);
  cx = (p(nodes1,1) + p(nodes2,1) + p(nodes3,1)) / 3;
  cy = (p(nodes1,2) + p(nodes2,2) + p(nodes3,2)) / 3;

  % Compute subdomain indices (ix, iy) in 1..nSubX / 1..nSubY.
  % Clamp to handle potential rounding at x=1 or y=1.
  ix = floor(cx * nSubX) + 1;
  iy = floor(cy * nSubY) + 1;
  ix = max(1, min(nSubX, ix));
  iy = max(1, min(nSubY, iy));

  elem_sub_id = ix + (iy-1)*nSubX; %flatten to column vector of the same dimension as <t>

  % -----------------------------
  % 4) Build per-subdomain element and node sets; build loc2glob
  % -----------------------------
  sub = repmat(struct(), ddm.Nsub, 1);

  for s = 1:ddm.Nsub

    % Coordinates of a subdomain in "subdomain" mesh. 
    % We define id = ix + (iy-1)*nSubX. 
    ix_s = mod(s-1, nSubX) + 1;
    iy_s = floor((s-1)/nSubX) + 1;

    %Subdomain's corners (x,y coordinates)
    x0 = (ix_s-1)/nSubX; x1 = ix_s/nSubX;
    y0 = (iy_s-1)/nSubY; y1 = iy_s/nSubY;

    elems = find(elem_sub_id == s); %returns indexes of elements, that lays ina subdomain <s>
    if isempty(elems)
      error('build_subdomains_structured: subdomain %d has no elements (unexpected for aligned partition).', s);
    end

    %extract all nodes (unique) that lays in a subdomain <s>
    nodes = unique(t(elems,:));
    nodes = nodes(:);

    % Convert nodes to global free DOFs (remove Dirichlet nodes by node2dof==0).
    % Mapping from local indexing to global one
    dofs = node2dof(nodes);
    dofs = dofs(dofs > 0);
    dofs = unique(dofs);
    dofs = sort(dofs);

    nloc = numel(dofs);

    % Map global free dof -> local index (0 if not present).
    glob2loc = zeros(ddm.nFree, 1);
    if nloc > 0
      glob2loc(dofs) = (1:nloc).';
    end

    %pack informations about one subdomain
    sub(s).id = s;
    sub(s).ix = ix_s;
    sub(s).iy = iy_s;
    sub(s).bbox = [x0 x1 y0 y1];
    sub(s).elems = elems;
    sub(s).nodes = nodes;
    sub(s).loc2glob = dofs;
    sub(s).glob2loc = glob2loc;
    sub(s).nloc = nloc;
  end
end
