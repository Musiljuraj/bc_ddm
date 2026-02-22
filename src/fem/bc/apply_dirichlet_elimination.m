function [Kff, Ff, free_nodes] = apply_dirichlet_elimination(K, F, dirichlet_nodes)
%APPLY_DIRICHLET_ELIMINATION  Enforce homogeneous Dirichlet BCs by restriction.
%
% Link to thesis:
%   Chapter 2.4.2–2.4.3 (Assembly and restriction of the discrete system).
%   Corresponds to the block decomposition (2.24)–(2.26) and the reduced system
%   formulated just on free (non-Dirichlet) degrees of freedom.
%
% Theoretical role:
% In the weak formulation of the Poisson problem with Dirichlet boundary
% conditions, the solution and test space is restricted to
%
%   V = { v in H^1(Ω) | v = 0 on Γ_D }.
%
% After discretization with P1 finite elements, this means that the degrees
% of freedom associated with Dirichlet nodes are fixed (here: u = 0) and
% therefore must be removed from the discrete system.
% Reordering the system conceptually gives:
%
%   [ K_FF  K_FD ] [ u_F ] = [ F_F ]
%   [ K_DF  K_DD ] [ u_D ]   [ F_D ]
%
% For homogeneous Dirichlet conditions, u_D = 0, and the reduced system (after removing Dirichlet nodes) is:
%
%   K_FF * u_F = F_F.
%
% Inputs:
%   K              : N×N global stiffness matrix (assembled on all nodes, no BCs applied yet)
%   F              : N×1 global load vector (assembled on all nodes, no BCs applied yet)
%   dirichlet_nodes: vector of global Dirichlet boundary nodes.
% Outputs:
%   Kff, Ff        : reduced matrix/vector on free nodes only.
%   free_nodes     : index list  of global indices of free degrees of freedom.
%
% Implementation outline:
%   1) Build the index set of free nodes as the complement of dirichlet_nodes.
%   2) Extract Kff = K(free,free) and Ff = F(free).
%   3) Return reduced system components for the subsequent linear solve.

% SOFTWARE DESIGN NOTES
% ---------------------
% - This function does NOT modify K or F in-place.
% - It performs only index manipulation and submatrix extraction.
% - This separation mirrors the theoretical idea that Dirichlet BCs restrict
%   the admissible function space, not the bilinear form itself.

  % ------------------------------------------------------------
  % 1) Input validation
  % ------------------------------------------------------------
  if nargin ~= 3
    error('apply_dirichlet_elimination: expected inputs (K, F, dirichlet_nodes).');
  end

  if ~ismatrix(K) || size(K,1) ~= size(K,2)
    error('apply_dirichlet_elimination: K must be a square matrix.');
  end

  if size(F,1) ~= size(K,1) || size(F,2) ~= 1
    error('apply_dirichlet_elimination: F must be a column vector of length size(K,1).');
  end

  % Type checks to avoid silent char/logical misuse.                   % ADDED
  if ~isnumeric(K) || ~isreal(K)                                       % ADDED
    error('apply_dirichlet_elimination: K must be a real numeric matrix.'); % ADDED
  end                                                                  % ADDED
  if ~isnumeric(F) || ~isreal(F)                                       % ADDED
    error('apply_dirichlet_elimination: F must be a real numeric column vector.'); % ADDED
  end                                                                  % ADDED

  N = size(K,1);

  % Ensure dirichlet_nodes is a vector
  dirichlet_nodes = dirichlet_nodes(:);

  % Harden index input: must be numeric, real, finite.                 % FIXED
  % (Prevents NaN/Inf being silently ignored by setdiff, and rejects    % ADDED
  %  char arrays like '23' that could otherwise be interpreted.)        % ADDED
  if islogical(dirichlet_nodes)                                        % ADDED
    error('apply_dirichlet_elimination: dirichlet_nodes must be an index vector (use find(mask)).'); % ADDED
  end                                                                  % ADDED
  if ~isnumeric(dirichlet_nodes) || ~isreal(dirichlet_nodes) || any(~isfinite(dirichlet_nodes)) % FIXED
    error('apply_dirichlet_elimination: dirichlet_nodes must be a real, finite numeric index vector.'); % FIXED
  end

  % Check index validity: integer-valued and within 1..N
  if any(dirichlet_nodes < 1) || any(dirichlet_nodes > N) || ...
     any(abs(dirichlet_nodes - round(dirichlet_nodes)) > 0)
    error('apply_dirichlet_elimination: dirichlet_nodes must be integer indices in 1..N.');
  end

  % Remove duplicates and sort (for robustness and reproducibility)
  dirichlet_nodes = unique(dirichlet_nodes);

  % ------------------------------------------------------------
  % 2) Determine free degrees of freedom
  % ------------------------------------------------------------
  % All global degrees of freedom:
  all_nodes = (1:N).';

  % Free nodes are those not subject to Dirichlet conditions
  free_nodes = setdiff(all_nodes, dirichlet_nodes);

  % Sanity check: there must be at least one free DOF
  if isempty(free_nodes)
    error('apply_dirichlet_elimination: no free degrees of freedom remain.');
  end

  % ------------------------------------------------------------
  % 3) Extract reduced system
  % ------------------------------------------------------------
  % This is the algebraic realization of the restriction:
  %
  %   Kff = K(F,F),  Ff = F(F),
  %
  % where F is the index set of free nodes.
  %
  % The resulting matrix Kff is symmetric positive definite
  % (for elliptic problems with Dirichlet BCs).
  %
  Kff = K(free_nodes, free_nodes);
  Ff  = F(free_nodes);

end