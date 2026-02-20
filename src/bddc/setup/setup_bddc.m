% ============================================================
% File: src/bddc/setup/setup_bddc.m
% ============================================================
function data = setup_bddc(data, opts)
%SETUP_BDDC  Precompute BDDC objects for the hat-space interface system.
%
% Baseline BDDC (Chapter 4.4):
%   - iteration space: assembled/hat interface space (size nHat)
%   - operator: A = R^T S R, where S = diag(S^(i)) in product space
%   - preconditioner: M^{-1} = R_D^T (T_sub + T_0) R_D
%       R_D = D R with multiplicity scaling D
%       T_sub = R_\Delta^T S_{\Delta\Delta}^{-1} R_\Delta
%       T_0   = \Psi (\Psi^T S \Psi)^{-1} \Psi^T
%
% Input:
%   data : problem data from build_problem_data(...)
%          must contain: sub, prod, primal
%   opts : optional struct
%          opts.tol_chol : SPD tolerance for chol (default 0)
%
% Output:
%   data : same struct, extended with fields in data.bddc

  if nargin < 2
    opts = struct();
  end
  if ~isfield(opts, 'tol_chol')
    opts.tol_chol = 0;
  end

  if ~isstruct(data) || ~isfield(data,'sub') || ~isfield(data,'prod') || ~isfield(data,'primal')
    error('setup_bddc: data must contain fields sub, prod, primal.');
  end

  sub = data.sub;
  prod = data.prod;
  primal = data.primal;
  nSub = numel(sub);

  % --- Build R (hat -> product) and multiplicity scaling weights ---
  R = build_assembly_operator_R(prod);  % sparse(nProd,nHat)
  omega = multiplicity_scaling(prod);   % (nProd x 1), omega=1/k on each hat group

  % --- Assemble product RHS g (condensed interface RHS) ---
  g_prod = zeros(prod.nProd, 1);
  for i = 1:nSub
    idx = sub(i).prod_idx(:);
    if isempty(idx)
      continue;
    end
    if ~isfield(sub(i),'g')
      error('setup_bddc: sub(i) missing field g (run setup_local_schur first).');
    end
    gi = sub(i).g(:);
    if numel(gi) ~= numel(idx)
      error('setup_bddc: sub(%d).g length mismatch with prod_idx.', i);
    end
    g_prod(idx) = gi;
  end

  % --- Local Schur block split (c,Delta) and Sdd factorisations ---
  Sdd_R = cell(nSub,1);
  for i = 1:nSub
    if ~isfield(sub(i),'S')
      error('setup_bddc: sub(i).S missing; build_problem_data uses opts.assemble_S=true.');
    end
    idx_d = sub(i).idx_d(:);
    if isempty(idx_d)
      Sdd_R{i} = [];
      continue;
    end
    Sdd = sub(i).S(idx_d, idx_d);
    [Rdd, pflag] = chol(Sdd);
    if pflag ~= 0
      error('setup_bddc: Sdd not SPD on subdomain %d (chol failed).', i);
    end
    Sdd_R{i} = Rdd; % upper triangular
  end

  % --- Build coarse basis Psi in product space ---
  nC = primal.nC;
  Psi = sparse(prod.nProd, nC);

  if nC > 0
    I = [];
    J = [];
    V = [];

    for j = 1:nC
      % Build column j by subdomain-wise energy-minimizing extension.
      col = zeros(prod.nProd, 1);

      for i = 1:nSub
        idx = sub(i).prod_idx(:);
        if isempty(idx)
          continue;
        end

        idx_c = sub(i).idx_c(:);
        idx_d = sub(i).idx_d(:);
        c_ids = sub(i).c_ids(:);

        nGi = numel(idx);
        if nGi == 0
          continue;
        end

        wloc = zeros(nGi, 1);

        % Set primal values (wc) on this subdomain for coarse id j.
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
          Sdc = sub(i).S(idx_d, idx_c);
          rhs = Sdc * wloc(idx_c);
          Rdd = Sdd_R{i};
          wD = -(Rdd \ (Rdd' \ rhs));
          wloc(idx_d) = wD;
        end

        col(idx) = wloc;
      end

      %transorm psi_j into sparse format
      nz = find(col ~= 0);
      if ~isempty(nz)
        I = [I; nz];
        J = [J; j*ones(numel(nz),1)];
        V = [V; col(nz)];
      end
    end

    Psi = sparse(I, J, V, prod.nProd, nC);
  end

  % --- Coarse matrix K0 = Psi^T S Psi and its factorization ---
  if nC > 0
    K0 = zeros(nC, nC);
    for j = 1:nC
      v = apply_blockdiag_S(sub, Psi(:,j));
      K0(:,j) = full(Psi' * v);
    end
    % Symmetrize (numerical) and factorize.
    K0 = 0.5 * (K0 + K0');
    [R0, pflag] = chol(K0);
    if pflag ~= 0
      error('setup_bddc: coarse matrix K0 is not SPD (chol failed).');
    end
  else
    K0 = zeros(0,0);
    R0 = [];
  end

  % --- Store ---
  b = struct();
  b.R = R;
  b.omega = omega;
  b.g_prod = g_prod;
  b.Sdd_R = Sdd_R;
  b.Psi = Psi;
  b.K0 = K0;
  b.K0_R = R0;

  data.bddc = b;
end