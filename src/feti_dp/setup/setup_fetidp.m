% ============================================================
% File: src/feti_dp/setup/setup_fetidp.m
% ============================================================
function data = setup_fetidp(data, opts)
%SETUP_FETIDP  Build all FETI-DP data required for multiplier-space PCG.
%
% Input:
%   data : struct from build_problem_data (Ch.2–3 pipeline)
%   opts : optional struct with fields:
%          - tol_rank (default 1e-12): rank threshold for Bd row reduction
%
% Output:
%   data : same struct extended with FETI-DP-specific fields:
%          - delta_idx_prod, nDeltaProd
%          - Scc/Scd/Sdc/Sdd (cell arrays), chol factors for inv(Sdd)
%          - Bd (full-row-rank reduced jump operator on Delta space), BdT
%          - weights D on Delta: DeltaWeights
%          - Kcc (coarse operator) + chol factor
%          - hc (coarse RHS part) for reconstruction
%          - delta_range per subdomain (mapping into packed Delta vector)

  if nargin < 2
    opts = struct();
  end
  if ~isfield(opts, 'tol_rank'); opts.tol_rank = 1e-12; end

  sub = data.sub;
  primal = data.primal;
  prod = data.prod;

  nSub = numel(sub);

  % ------------------------------------------------------------
  % Delta indices in product space
  % ------------------------------------------------------------
  delta_idx = primal.delta_idx(:);           % product indices of W_Delta
  data.delta_idx_prod = delta_idx;
  data.nDeltaProd = numel(delta_idx);

  % Build per-subdomain ranges into the packed delta vector.
  % For each subdomain i:
  %   sub(i).prod_idx_d is a list of product indices belonging to that subdomain's Delta subset.
  % We map these product indices to their position inside delta_idx (packed order).
  prod2delta = primal.prod2delta;            % size nProd, 0 if primal, else 1..nDeltaProd
  data.delta_range = cell(nSub,1);
  for i = 1:nSub
    if ~isfield(sub(i), 'prod_idx_d')
      error('setup_fetidp: missing sub(i).prod_idx_d (run build_primal_maps).');
    end
    if isempty(sub(i).prod_idx_d)
      data.delta_range{i} = zeros(0,1);
    else
      data.delta_range{i} = prod2delta(sub(i).prod_idx_d(:));
    end
  end

  % ------------------------------------------------------------
  % Local Schur blocks split into [c, Delta] in LOCAL gamma ordering
  % ------------------------------------------------------------
  data.Scc = cell(nSub,1);
  data.Scd = cell(nSub,1);
  data.Sdc = cell(nSub,1);
  data.Sdd = cell(nSub,1);
  data.Sdd_R = cell(nSub,1);

  for i = 1:nSub
    if ~isfield(sub(i), 'S') || isempty(sub(i).S)
      error('setup_fetidp: sub(i).S missing. Run setup_local_schur with opts.assemble_S=true.');
    end
    S = sub(i).S;

    idx_c = sub(i).idx_c(:);
    idx_d = sub(i).idx_d(:);

    % Handle empty interface or empty delta safely.
    if isempty(idx_c)
      data.Scc{i} = sparse(0,0);
      data.Scd{i} = sparse(0, numel(idx_d));
      data.Sdc{i} = sparse(numel(idx_d), 0);
    else
      data.Scc{i} = S(idx_c, idx_c);
      data.Scd{i} = S(idx_c, idx_d);
      data.Sdc{i} = S(idx_d, idx_c);
    end

    if isempty(idx_d)
      data.Sdd{i} = sparse(0,0);
      data.Sdd_R{i} = [];
    else
      data.Sdd{i} = S(idx_d, idx_d);
      [R,pflag] = chol(full(data.Sdd{i})); 
      if pflag ~= 0
        error('setup_fetidp: Sdd is not SPD on subdomain %d (pflag=%d).', i, pflag);
      end
      data.Sdd_R{i} = R; % upper
    end
  end

  % ------------------------------------------------------------
  % Build reduced jump operator Bd on Delta space, then row-reduce to full row rank
  % ------------------------------------------------------------
  B = build_jump_operator_B(prod);
  Bd_full = B(:, delta_idx);                 % restrict columns to Delta

  % Remove dependent rows via QR on Bd_full'
  % In Octave: qr(A,0) with pivoting returns E as permutation indices.
  if isempty(Bd_full)
    Bd = sparse(0, numel(delta_idx));
    ind_rows = zeros(0,1);
  else
    %[~, Rq, E] = qr(Bd_full', 0);
    %d = abs(diag(Rq));
    %r = sum(d > opts.tol_rank);
    %ind_rows = E(1:r);
    %Bd = Bd_full(ind_rows, :);
    [~, Rq, E] = qr(Bd_full', 0);
    d = abs(diag(Rq));
    r = sum(d > opts.tol_rank);

    if isvector(E)
      p = E(:);
    else
      % E is a permutation matrix; convert to permutation vector p with A(:,p) = Q*R
      [~, p] = max(E, [], 1);
      p = p(:);
    end

    ind_rows = p(1:r);
    Bd = Bd_full(ind_rows, :);
  end

  data.Bd = sparse(Bd);
  data.BdT = data.Bd';
  data.lambda_rows = ind_rows(:);
  data.nLambda = size(data.Bd,1);

  % ------------------------------------------------------------
  % Multiplicity weights on product space, restricted to Delta
  % ------------------------------------------------------------
  omega_prod = multiplicity_scaling(prod);
  data.DeltaWeights = omega_prod(delta_idx);

  % ------------------------------------------------------------
  % Build coarse operator Kcc* and factorization (assembled; small)
  % Kcc = sum_i A_i^T (Scc - Scd*Sdd^{-1}*Sdc) A_i
  % Here A_i is the boolean injection to global coarse ids via sub(i).c_ids.
  % Also build hc = sum_i A_i^T (gc - Scd*Sdd^{-1}*gd)
  % ------------------------------------------------------------
  nC = primal.nC; %Number of global coarse/primal dofs.
  Kcc = sparse(nC, nC); %Initialize sparse global coarse matrix.
  hc = zeros(nC,1); %Initialize coarse RHS vector.

  for i = 1:nSub
    c_ids = sub(i).c_ids(:); %Local-to-global coarse id map for this subdomain’s primal dofs.
    if isempty(c_ids)
      continue;
    end

    Scc = data.Scc{i};
    Scd = data.Scd{i};
    Sdc = data.Sdc{i};

    % Local RHS split from sub(i).g (condensed interface RHS in local gamma ordering)
    if ~isfield(sub(i), 'g')
      error('setup_fetidp: sub(i).g missing. Ensure setup_local_schur ran.');
    end
    gi = sub(i).g(:); %Local condensed RHS vector.
    gc_i = gi(sub(i).idx_c(:)); %Extract primal component
    gd_i = gi(sub(i).idx_d(:)); %extract delta component

    if isempty(sub(i).idx_d)
      localK = Scc;
      localh = gc_i;
    else
      R = data.Sdd_R{i};  %get already computed Cholesky factor of Sdd
      % X = Sdd^{-1} * Sdc
      X = R \ (R' \ full(Sdc));
      localK = Scc - Scd * X; %compute local Kcc^i

      % y = Sdd^{-1} * gd
      y = R \ (R' \ full(gd_i));
      localh = gc_i - Scd * y;%compute local h_c^i
    end

    % Assemble into global coarse system
    Kcc(c_ids, c_ids) = Kcc(c_ids, c_ids) + sparse(localK); 
    hc(c_ids) = hc(c_ids) + localh;
  end

  data.Kcc = Kcc;
  data.hc = hc;

  if nC > 0
    [Rcc,pflag] = chol(full(Kcc)); %compute Cholesky factor of K_cc, once
    if pflag ~= 0
      error('setup_fetidp: Kcc is not SPD (pflag=%d).', pflag);
    end
    data.Kcc_R = Rcc; % upper
  else
    data.Kcc_R = [];
  end

  % Keep the updated sub array in data
  data.sub = sub;
end