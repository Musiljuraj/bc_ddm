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

  % ----------------------------
  % Input validation / defaults
  % ----------------------------
  % ADDED: strict nargin handling and opts validation
  if nargin < 2 || isempty(opts)
    opts = struct();
  end
  if ~isstruct(opts)
    error('setup_fetidp: opts must be a struct when provided.');
  end
  if ~isfield(opts, 'tol_rank'); opts.tol_rank = 1e-12; end
  if ~(isnumeric(opts.tol_rank) && isreal(opts.tol_rank) && isscalar(opts.tol_rank) && isfinite(opts.tol_rank) && opts.tol_rank >= 0)
    error('setup_fetidp: opts.tol_rank must be a real, finite, nonnegative scalar.');
  end

  % ADDED: data contract checks
  if ~isstruct(data)
    error('setup_fetidp: data must be a struct.');
  end
  if ~isfield(data,'sub') || ~isfield(data,'primal') || ~isfield(data,'prod')
    error('setup_fetidp: data must contain fields sub, primal, prod.');
  end
  if ~isstruct(data.sub)
    error('setup_fetidp: data.sub must be a struct array.');
  end
  if ~isstruct(data.primal)
    error('setup_fetidp: data.primal must be a struct.');
  end

  sub = data.sub;
  primal = data.primal;
  prod = data.prod;

  nSub = numel(sub);
  if nSub < 1
    error('setup_fetidp: data.sub must contain at least one subdomain.');
  end

  % ADDED: validate primal fields used by this routine
  if ~isfield(primal,'delta_idx') || ~isfield(primal,'prod2delta') || ~isfield(primal,'nC')
    error('setup_fetidp: data.primal must contain delta_idx, prod2delta, nC.');
  end
  if ~(isnumeric(primal.nC) && isreal(primal.nC) && isscalar(primal.nC) && isfinite(primal.nC) && primal.nC >= 0 && abs(primal.nC - round(primal.nC)) == 0)
    error('setup_fetidp: primal.nC must be a nonnegative integer scalar.');
  end

  prod2delta = primal.prod2delta;
  if ~(isnumeric(prod2delta) && isreal(prod2delta))
    error('setup_fetidp: primal.prod2delta must be real numeric.');
  end
  if any(~isfinite(prod2delta(:)))
    error('setup_fetidp: primal.prod2delta must be finite.');
  end
  if any(abs(prod2delta(:) - round(prod2delta(:))) > 0)
    error('setup_fetidp: primal.prod2delta must be integer-valued.');
  end
  nProd = numel(prod2delta);

  delta_idx = primal.delta_idx(:);
  if ~(isnumeric(delta_idx) && isreal(delta_idx))
    error('setup_fetidp: primal.delta_idx must be real numeric indices.');
  end
  if any(~isfinite(delta_idx(:)))
    error('setup_fetidp: primal.delta_idx must be finite.');
  end
  if any(abs(delta_idx(:) - round(delta_idx(:))) > 0)
    error('setup_fetidp: primal.delta_idx must be integer-valued.');
  end
  if any(delta_idx < 1) || any(delta_idx > nProd)
    error('setup_fetidp: primal.delta_idx entries out of range for primal.prod2delta.');
  end

  % ------------------------------------------------------------
  % Delta indices in product space
  % ------------------------------------------------------------
  % CHANGED: moved delta_idx definition above (validation needs nProd)
  data.delta_idx_prod = delta_idx;
  data.nDeltaProd = numel(delta_idx);

  % ADDED: consistency check between delta_idx (packed order) and prod2delta mapping
  if data.nDeltaProd > 0
    mapped = prod2delta(delta_idx);
    if any(mapped < 1) || any(mapped > data.nDeltaProd)
      error('setup_fetidp: primal.prod2delta(delta_idx) must lie in 1..nDeltaProd.');
    end
    if ~isequal(mapped(:), (1:data.nDeltaProd).')
      error('setup_fetidp: inconsistent primal.delta_idx / primal.prod2delta mapping (expected prod2delta(delta_idx(k)) == k).');
    end
  end

  % Build per-subdomain ranges into the packed delta vector.
  data.delta_range = cell(nSub,1);
  for i = 1:nSub
    if ~isfield(sub(i), 'prod_idx_d')
      error('setup_fetidp: missing sub(i).prod_idx_d (run build_primal_maps).');
    end
    pid = sub(i).prod_idx_d(:);

    % ADDED: validate prod_idx_d indices before mapping
    if ~(isnumeric(pid) && isreal(pid))
      error('setup_fetidp: sub(%d).prod_idx_d must be real numeric indices.', i);
    end
    if any(~isfinite(pid))
      error('setup_fetidp: sub(%d).prod_idx_d must be finite.', i);
    end
    if any(abs(pid - round(pid)) > 0)
      error('setup_fetidp: sub(%d).prod_idx_d must be integer-valued.', i);
    end
    if any(pid < 1) || any(pid > nProd)
      error('setup_fetidp: sub(%d).prod_idx_d entries out of range for primal.prod2delta.', i);
    end

    if isempty(pid)
      data.delta_range{i} = zeros(0,1);
    else
      dr = prod2delta(pid);

      % ADDED: disallow primal indices inside prod_idx_d (mapping to 0)
      if any(dr == 0)
        error('setup_fetidp: sub(%d).prod_idx_d contains non-Delta (primal) product indices (prod2delta==0).', i);
      end
      if any(dr < 1) || any(dr > data.nDeltaProd)
        error('setup_fetidp: sub(%d).prod_idx_d maps outside 1..nDeltaProd.', i);
      end
      data.delta_range{i} = dr;
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

    % ADDED: validate S is square and finite real
    if ~(isnumeric(S) && isreal(S))
      error('setup_fetidp: sub(%d).S must be a real numeric matrix.', i);
    end
    if any(~isfinite(S(:)))
      error('setup_fetidp: sub(%d).S must be finite.', i);
    end
    if size(S,1) ~= size(S,2)
      error('setup_fetidp: sub(%d).S must be square.', i);
    end
    nGamma = size(S,1);

    % ADDED: validate idx_c/idx_d exist and are indices
    if ~isfield(sub(i),'idx_c') || ~isfield(sub(i),'idx_d')
      error('setup_fetidp: sub(%d) missing idx_c/idx_d (run build_primal_maps / interface split).', i);
    end
    idx_c = sub(i).idx_c(:);
    idx_d = sub(i).idx_d(:);

    if ~(isnumeric(idx_c) && isreal(idx_c) && isnumeric(idx_d) && isreal(idx_d))
      error('setup_fetidp: sub(%d).idx_c/idx_d must be real numeric indices.', i);
    end
    if any(~isfinite([idx_c; idx_d]))
      error('setup_fetidp: sub(%d).idx_c/idx_d must be finite.', i);
    end
    if any(abs([idx_c; idx_d] - round([idx_c; idx_d])) > 0)
      error('setup_fetidp: sub(%d).idx_c/idx_d must be integer-valued.', i);
    end
    if any(idx_c < 1) || any(idx_c > nGamma) || any(idx_d < 1) || any(idx_d > nGamma)
      error('setup_fetidp: sub(%d).idx_c/idx_d entries out of range for S.', i);
    end
    % ADDED: enforce disjointness of primal and delta local index sets
    if ~isempty(intersect(idx_c, idx_d))
      error('setup_fetidp: sub(%d).idx_c and idx_d must be disjoint.', i);
    end

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

      % CHANGED: keep original chol(full()) approach, but validate failure clearly
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
  % ADDED: dependency checks (fail fast with a clear message)
  if exist('build_jump_operator_B','file') ~= 2
    error('setup_fetidp: build_jump_operator_B not found on path.');
  end
  if exist('multiplicity_scaling','file') ~= 2
    error('setup_fetidp: multiplicity_scaling not found on path.');
  end

  B = build_jump_operator_B(prod);
  Bd_full = B(:, delta_idx);                 % restrict columns to Delta

  if isempty(Bd_full)
    Bd = sparse(0, numel(delta_idx));
    ind_rows = zeros(0,1);
  else
    [~, Rq, E] = qr(Bd_full', 0);
    d = abs(diag(Rq));
    r = sum(d > opts.tol_rank);

    % ADDED: guard r bounds (defensive)
    r = max(0, min(r, size(Bd_full,1)));

    if isvector(E)
      p = E(:);
    else
      [~, p] = max(E, [], 1);
      p = p(:);
    end

    if r > numel(p)
      error('setup_fetidp: internal pivot selection exceeded available rows (r=%d, n=%d).', r, numel(p));
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

  % ADDED: validate omega_prod indexing and positivity
  if ~(isnumeric(omega_prod) && isreal(omega_prod)) || any(~isfinite(omega_prod(:)))
    error('setup_fetidp: multiplicity_scaling(prod) must return a finite real numeric vector.');
  end
  if numel(omega_prod) < nProd
    error('setup_fetidp: multiplicity_scaling(prod) returned vector shorter than nProd.');
  end
  data.DeltaWeights = omega_prod(delta_idx);
  if any(~isfinite(data.DeltaWeights(:))) || any(data.DeltaWeights(:) <= 0)
    error('setup_fetidp: DeltaWeights must be positive and finite.');
  end

  % ------------------------------------------------------------
  % Build coarse operator Kcc* and factorization (assembled; small)
  % ------------------------------------------------------------
  nC = primal.nC; %Number of global coarse/primal dofs.
  Kcc = sparse(nC, nC);
  hc = zeros(nC,1);

  for i = 1:nSub
    % ADDED: validate c_ids exists if idx_c is used (assembly correctness)
    if ~isfield(sub(i),'c_ids')
      error('setup_fetidp: sub(%d).c_ids missing (run build_primal_maps).', i);
    end
    c_ids = sub(i).c_ids(:);

    if ~(isnumeric(c_ids) && isreal(c_ids))
      error('setup_fetidp: sub(%d).c_ids must be real numeric indices.', i);
    end
    if any(~isfinite(c_ids))
      error('setup_fetidp: sub(%d).c_ids must be finite.', i);
    end
    if any(abs(c_ids - round(c_ids)) > 0)
      error('setup_fetidp: sub(%d).c_ids must be integer-valued.', i);
    end
    if any(c_ids < 1) || any(c_ids > nC)
      error('setup_fetidp: sub(%d).c_ids out of range 1..nC.', i);
    end

    % ADDED: require c_ids length consistent with idx_c (common contract)
    if numel(c_ids) ~= numel(sub(i).idx_c(:))
      error('setup_fetidp: sub(%d) mismatch: numel(c_ids) must equal numel(idx_c).', i);
    end

    if isempty(c_ids)
      continue;
    end

    Scc = data.Scc{i};
    Scd = data.Scd{i};
    Sdc = data.Sdc{i};

    if ~isfield(sub(i), 'g')
      error('setup_fetidp: sub(i).g missing. Ensure setup_local_schur ran.');
    end
    gi = sub(i).g(:);

    % ADDED: validate g length consistent with S (local gamma ordering)
    if ~(isnumeric(gi) && isreal(gi)) || any(~isfinite(gi))
      error('setup_fetidp: sub(%d).g must be finite real numeric.', i);
    end
    % size(S,1) was validated earlier per-subdomain; recover via blocks
    nGamma_i = numel(sub(i).idx_c(:)) + numel(sub(i).idx_d(:));
    if numel(gi) < max([sub(i).idx_c(:); sub(i).idx_d(:); 0])
      error('setup_fetidp: sub(%d).g length too small for idx_c/idx_d.', i);
    end

    gc_i = gi(sub(i).idx_c(:));
    gd_i = gi(sub(i).idx_d(:));

    if isempty(sub(i).idx_d)
      localK = Scc;
      localh = gc_i;
    else
      R = data.Sdd_R{i};

      X = R \ (R' \ full(Sdc));
      localK = Scc - Scd * X;

      y = R \ (R' \ full(gd_i));
      localh = gc_i - Scd * y;
    end

    Kcc(c_ids, c_ids) = Kcc(c_ids, c_ids) + sparse(localK);
    hc(c_ids) = hc(c_ids) + localh;
  end

  data.Kcc = Kcc;
  data.hc = hc;

  if nC > 0
    [Rcc,pflag] = chol(full(Kcc));
    if pflag ~= 0
      error('setup_fetidp: Kcc is not SPD (pflag=%d).', pflag);
    end
    data.Kcc_R = Rcc;
  else
    data.Kcc_R = [];
  end

  data.sub = sub;
end