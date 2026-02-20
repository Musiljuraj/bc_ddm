function sub = setup_local_schur(sub, opts)
%SETUP_LOCAL_SCHUR  Prepare local interior solves and reduced RHS for Schur complement formulation.
%
% Link to thesis:
%   Chapter 3.2.2–3.2.4 (interior elimination, Schur complement),
%   equations (3.19)–(3.21), (3.29), (3.33)–(3.34).
%
% For each subdomain i, given the block system
%   [K_II  K_Ig] [u_I] = [f_I]
%   [K_gI  K_gg] [u_g]   [f_g]
% we precompute:
%   - a Cholesky factorization of K_II (local interior solver),
%   - the reduced RHS:
%       g = f_g - K_gI * (K_II^{-1} * f_I)
%
% Optionally, an explicit Schur complement matrix can be assembled:
%       S = K_gg - K_gI * (K_II^{-1} * K_Ig)
% but for large problems matrix-free application is preferred.
%
% Inputs:
%   sub  : subdomain struct array with fields K_II, K_Ig, K_gI, K_gg, f_I, f_g.
%   opts : optional struct with fields:
%          .assemble_S (default false) : if true, compute and store explicit S.
%
% Output:
%   sub : updated with fields:
%       .R_II : upper Cholesky factor of K_II (empty if nI=0), useh heavily later for computation with Schur complements
%       .g    : reduced RHS on interface DOFs
%       .S    : explicit Schur complement (only if opts.assemble_S=true)

  if nargin < 1 || nargin > 2
    error('setup_local_schur: expected inputs (sub[,opts]).');
  end
  if nargin == 1
    opts = struct();
  end
  if ~isfield(opts,'assemble_S')
    opts.assemble_S = false;
  end

  nSub = numel(sub);

  for i = 1:nSub
    if ~isfield(sub(i),'K_II')
      error('setup_local_schur: sub(%d) missing block matrices (run extract_subdomain_blocks first).', i);
    end

    nI = size(sub(i).K_II,1);
    nG = size(sub(i).K_gg,1);

    % Ensure RHS fields exist, even if empty.
    if ~isfield(sub(i),'f_I') || ~isfield(sub(i),'f_g')
      error('setup_local_schur: sub(%d) missing RHS blocks (run extract_subdomain_blocks first).', i);
    end

    if nI == 0
      sub(i).R_II = [];
      sub(i).g = sub(i).f_g;

      if opts.assemble_S
        sub(i).S = sub(i).K_gg;
      end
      continue;
    end

    % Cholesky factorization (SPD).
    % chol returns upper triangular R such that R'*R = K_II.
    [R, pflag] = chol(sub(i).K_II);
    if pflag ~= 0
      error('setup_local_schur: chol failed for sub(%d).K_II (matrix not SPD?).', i);
    end
    sub(i).R_II = R;

    % Reduced RHS: g = f_g - K_gI * (K_II^{-1} f_I)
    z = R \ (R' \ sub(i).f_I); %Cholesky solve
    sub(i).g = sub(i).f_g - sub(i).K_gI * z;

    if opts.assemble_S
      X = R \ (R' \ sub(i).K_Ig);
      sub(i).S = sub(i).K_gg - sub(i).K_gI * X;
    end
  end
end
