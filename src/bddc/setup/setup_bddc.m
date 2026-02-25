function data = setup_bddc(data, opts)
%SETUP_BDDC  Precompute BDDC objects for the hat-space interface system.
%
% Input:
%   data : problem data (struct) with fields: sub, prod, primal
%   opts : optional struct, supports:
%          opts.tol_chol : nonnegative diagonal shift used only if chol fails (default 0)
%
% Output:
%   data : same struct, extended with data.bddc

  %======================================================================
  %% HARDENING CHANGE: robust opts handling + validation (non-struct, missing fields, bad tol)
  %======================================================================
  if nargin < 2 || isempty(opts)
    opts = struct();
  end
  if ~isstruct(opts)
    error('setup_bddc:opts', 'opts must be a struct.');
  end
  if ~isfield(opts, 'tol_chol') || isempty(opts.tol_chol)
    opts.tol_chol = 0;
  end
  if ~(isnumeric(opts.tol_chol) && isscalar(opts.tol_chol) && isfinite(opts.tol_chol) && opts.tol_chol >= 0)
    error('setup_bddc:opts', 'opts.tol_chol must be a nonnegative finite scalar.');
  end

  %======================================================================
  %% HARDENING CHANGE: top-level input validation with explicit error IDs
  %======================================================================
  if ~isstruct(data) || ~isfield(data,'sub') || ~isfield(data,'prod') || ~isfield(data,'primal')
    error('setup_bddc:data', 'data must be a struct containing fields sub, prod, primal.');
  end

  sub    = data.sub;
  prod   = data.prod;
  primal = data.primal;

  if ~isstruct(sub)
    error('setup_bddc:data', 'data.sub must be a struct array.');
  end
  if ~isstruct(prod) || ~isfield(prod,'nProd') || ~isfield(prod,'nHat')
    error('setup_bddc:data', 'data.prod must contain fields nProd and nHat.');
  end
  if ~isstruct(primal) || ~isfield(primal,'nC')
    error('setup_bddc:data', 'data.primal must contain field nC.');
  end

  nSub  = numel(sub);
  nProd = prod.nProd;
  nHat  = prod.nHat;
  nC    = primal.nC;

  %======================================================================
  %% HARDENING CHANGE: enforce integer sizes for nProd/nHat/nC
  %======================================================================
  if ~(isnumeric(nProd) && isscalar(nProd) && nProd >= 0 && floor(nProd) == nProd)
    error('setup_bddc:data', 'prod.nProd must be a nonnegative integer scalar.');
  end
  if ~(isnumeric(nHat) && isscalar(nHat) && nHat >= 0 && floor(nHat) == nHat)
    error('setup_bddc:data', 'prod.nHat must be a nonnegative integer scalar.');
  end
  if ~(isnumeric(nC) && isscalar(nC) && nC >= 0 && floor(nC) == nC)
    error('setup_bddc:data', 'primal.nC must be a nonnegative integer scalar.');
  end

  %======================================================================
  %% HARDENING CHANGE: validate dependency outputs (R/omega sizes, omega positivity/finite)
  %======================================================================
  R = build_assembly_operator_R(prod);  % expected sparse(nProd,nHat)
  omega = multiplicity_scaling(prod);   % expected (nProd x 1)

  if ~(isnumeric(R) && all(size(R) == [nProd, nHat]))
    error('setup_bddc:R', 'build_assembly_operator_R(prod) must return size [nProd,nHat].');
  end
  if ~issparse(R)
    R = sparse(R);
  end
  if ~(isnumeric(omega) && isvector(omega) && numel(omega) == nProd)
    error('setup_bddc:omega', 'multiplicity_scaling(prod) must return a vector of length nProd.');
  end
  omega = omega(:);
  if any(~isfinite(omega)) || any(omega <= 0)
    error('setup_bddc:omega', 'omega must be positive and finite.');
  end

  %======================================================================
  %% Assemble product RHS g_prod
  %% HARDENING CHANGE: detect out-of-range indices + overlapping prod_idx blocks
  %======================================================================
  g_prod = zeros(nProd, 1);
  seen = false(max(nProd,1), 1); % overlap detection

  for i = 1:nSub
    if ~isfield(sub(i),'prod_idx')
      error('setup_bddc:data', 'sub(%d) missing field prod_idx.', i);
    end
    idx = sub(i).prod_idx(:);
    if isempty(idx)
      continue;
    end

    % HARDENING CHANGE: prod_idx integer + bounds check
    if ~(isnumeric(idx) && all(isfinite(idx)) && all(floor(idx) == idx))
      error('setup_bddc:data', 'sub(%d).prod_idx must contain integer indices.', i);
    end
    if any(idx < 1) || any(idx > nProd)
      error('setup_bddc:data', 'sub(%d).prod_idx out of range 1..prod.nProd.', i);
    end

    % HARDENING CHANGE: overlap check (product indices should be disjoint)
    if any(seen(idx))
      error('setup_bddc:data', 'Overlapping prod_idx detected (product space indices must be disjoint).');
    end
    seen(idx) = true;

    if ~isfield(sub(i),'g')
      error('setup_bddc:data', 'sub(%d) missing field g (run setup_local_schur first).', i);
    end
    gi = sub(i).g(:);

    % HARDENING CHANGE: length and finiteness check
    if ~(isnumeric(gi) && numel(gi) == numel(idx))
      error('setup_bddc:data', 'sub(%d).g length mismatch with prod_idx.', i);
    end
    if any(~isfinite(gi))
      error('setup_bddc:data', 'sub(%d).g contains NaN/Inf.', i);
    end

    g_prod(idx) = gi;
  end

  %======================================================================
  %% Factorize local Sdd blocks
  %% HARDENING CHANGE: validate S shape, idx_c/idx_d/c_ids ranges/disjointness,
  %%                   chol with p-flag + optional diagonal shift + informative error
  %======================================================================
  Sdd_R = cell(nSub,1);

  for i = 1:nSub
    idx = sub(i).prod_idx(:);
    nGi = numel(idx);
    if nGi == 0
      Sdd_R{i} = [];
      continue;
    end

    if ~isfield(sub(i),'S') || isempty(sub(i).S)
      error('setup_bddc:data', 'sub(%d).S missing/empty; ensure local Schur S is assembled.', i);
    end
    Si = sub(i).S;

    % HARDENING CHANGE: S must match prod_idx size
    if ~(isnumeric(Si) && size(Si,1) == nGi && size(Si,2) == nGi)
      error('setup_bddc:data', 'sub(%d).S must be size [numel(prod_idx), numel(prod_idx)].', i);
    end

    if ~isfield(sub(i),'idx_c') || ~isfield(sub(i),'idx_d') || ~isfield(sub(i),'c_ids')
      error('setup_bddc:data', 'sub(%d) must contain idx_c, idx_d, c_ids.', i);
    end

    idx_c = sub(i).idx_c(:);
    idx_d = sub(i).idx_d(:);
    c_ids = sub(i).c_ids(:);

    % HARDENING CHANGE: local index validation
    if ~(isnumeric(idx_c) && all(isfinite(idx_c)) && all(floor(idx_c) == idx_c))
      error('setup_bddc:data', 'sub(%d).idx_c must be integer indices.', i);
    end
    if ~(isnumeric(idx_d) && all(isfinite(idx_d)) && all(floor(idx_d) == idx_d))
      error('setup_bddc:data', 'sub(%d).idx_d must be integer indices.', i);
    end
    if any(idx_c < 1) || any(idx_c > nGi) || any(idx_d < 1) || any(idx_d > nGi)
      error('setup_bddc:data', 'sub(%d).idx_c/idx_d out of local range 1..numel(prod_idx).', i);
    end
    if ~isempty(intersect(idx_c, idx_d))
      error('setup_bddc:data', 'sub(%d): idx_c and idx_d must be disjoint.', i);
    end

    % HARDENING CHANGE: c_ids validation + range check when nC > 0
    if ~(isnumeric(c_ids) && all(isfinite(c_ids)) && all(floor(c_ids) == c_ids))
      error('setup_bddc:data', 'sub(%d).c_ids must be integer coarse ids.', i);
    end
    if ~isempty(idx_c) && numel(c_ids) ~= numel(idx_c)
      error('setup_bddc:data', 'sub(%d): numel(c_ids) must match numel(idx_c).', i);
    end
    if ~isempty(c_ids)
      if nC == 0
        error('setup_bddc:data', 'sub(%d): nonempty c_ids but primal.nC==0.', i);
      end
      if any(c_ids < 1) || any(c_ids > nC)
        error('setup_bddc:data', 'sub(%d).c_ids out of range 1..primal.nC.', i);
      end
    end

    if isempty(idx_d)
      Sdd_R{i} = [];
      continue;
    end

    Sdd = Si(idx_d, idx_d);

    % HARDENING CHANGE: symmetrize block (guards tiny asymmetry)
    Sdd = 0.5 * (Sdd + Sdd');

    % HARDENING CHANGE: chol with p-flag and contextual error
    [Rdd, pflag] = chol(Sdd);
    if pflag ~= 0 && opts.tol_chol > 0
      Sdd = Sdd + opts.tol_chol * speye(size(Sdd,1));
      [Rdd, pflag] = chol(Sdd);
    end
    if pflag ~= 0
      error('setup_bddc:chol', 'Sdd not SPD on subdomain %d (chol failed).', i);
    end

    Sdd_R{i} = Rdd; % upper triangular
  end

  %======================================================================
  %% Build coarse basis Psi in product space
  %% HARDENING CHANGE: handle nC==0 gracefully with correct empty shapes
  %======================================================================
  Psi = sparse(nProd, nC);

  if nC > 0
    I = [];
    J = [];
    V = [];

    for j = 1:nC
      col = zeros(nProd, 1);

      for i = 1:nSub
        idx = sub(i).prod_idx(:);
        if isempty(idx)
          continue;
        end

        idx_c = sub(i).idx_c(:);
        idx_d = sub(i).idx_d(:);
        c_ids = sub(i).c_ids(:);
        nGi = numel(idx);

        wloc = zeros(nGi, 1);

        % Set primal values on this subdomain for coarse id j.
        if ~isempty(idx_c)
          wc = zeros(numel(idx_c), 1);
          for k = 1:numel(c_ids)
            if c_ids(k) == j
              wc(k) = 1.0;
            end
          end
          wloc(idx_c) = wc;
        end

        % Energy-minimizing extension into Delta: wD = -Sdd^{-1} * Sdc * wc.
        if ~isempty(idx_d) && ~isempty(idx_c)
          Si = sub(i).S;
          Sdc = Si(idx_d, idx_c);
          rhs = Sdc * wloc(idx_c);

          Rdd = Sdd_R{i};
          wD = -(Rdd \ (Rdd' \ rhs)); % Sdd wD = -rhs
          wloc(idx_d) = wD;
        end

        col(idx) = wloc;
      end

      nz = find(col ~= 0);
      if ~isempty(nz)
        I = [I; nz];
        J = [J; j*ones(numel(nz),1)];
        V = [V; col(nz)];
      end
    end

    Psi = sparse(I, J, V, nProd, nC);
  end

  %======================================================================
  %% Coarse matrix K0 = Psi^T S Psi
  %% HARDENING CHANGE: validate apply_blockdiag_S output, symmetrize K0,
  %%                   chol with p-flag + optional shift + informative error
  %======================================================================
  if nC > 0
    if exist('apply_blockdiag_S','file') ~= 2
      error('setup_bddc:deps', 'apply_blockdiag_S not found on path.');
    end

    K0 = zeros(nC, nC);
    for j = 1:nC
      v = apply_blockdiag_S(sub, Psi(:,j));

      % HARDENING CHANGE: check size/type returned by apply_blockdiag_S
      if ~(isnumeric(v) && isvector(v) && numel(v) == nProd)
        error('setup_bddc:K0', 'apply_blockdiag_S returned invalid vector size.');
      end

      K0(:,j) = full(Psi' * v(:));
    end

    % HARDENING CHANGE: symmetrize K0 (numerical)
    K0 = 0.5 * (K0 + K0');

    [R0, pflag] = chol(K0);
    if pflag ~= 0 && opts.tol_chol > 0
      K0s = K0 + opts.tol_chol * eye(nC);
      [R0, pflag] = chol(K0s);
      if pflag == 0
        % HARDENING CHANGE: keep K0 consistent with actual factorization used
        K0 = K0s;
      end
    end
    if pflag ~= 0
      error('setup_bddc:chol', 'Coarse matrix K0 is not SPD (chol failed).');
    end
  else
    % HARDENING CHANGE: explicit empty outputs when nC==0
    K0 = zeros(0,0);
    R0 = [];
  end

  %======================================================================
  %% Store outputs
  %======================================================================
  b = struct();
  b.R      = R;
  b.omega  = omega;
  b.g_prod = g_prod;
  b.Sdd_R  = Sdd_R;
  b.Psi    = Psi;
  b.K0     = K0;
  b.K0_R   = R0;

  data.bddc = b;
end