% ============================================================
% File: src/feti_dp/operators/solve_tildeS.m
% ============================================================
function out = solve_tildeS(rc, rD, data)
%SOLVE_TILDES  Apply \widetilde{S}^{-1} to RHS split into primal and delta parts.
%
% Output "out" struct:
%   out.u_c        : coarse vector (nC)
%   out.wc_local{i}: local primal values on subdomain i (length nC_i)
%   out.wD         : packed Delta vector (nDeltaProd)

  % FIXED: require data argument (use [] for defaults of rc/rD)
  if nargin < 3
    error('solve_tildeS:argchk', ...
          'solve_tildeS requires rc, rD, and data (use [] for default rc/rD).');
  end

  % ADDED: minimal validation of data struct and required fields
  if ~isstruct(data)
    error('solve_tildeS:badData', 'data must be a struct.');
  end
  if ~isfield(data,'sub') || ~isstruct(data.sub)
    error('solve_tildeS:badData', 'data.sub must be a struct array.');
  end
  if ~isfield(data,'primal') || ~isstruct(data.primal) || ~isfield(data.primal,'nC')
    error('solve_tildeS:badData', 'data.primal.nC is required.');
  end
  if ~isfield(data,'nDeltaProd')
    error('solve_tildeS:badData', 'data.nDeltaProd is required.');
  end
  if ~isfield(data,'delta_range') || ~iscell(data.delta_range)
    error('solve_tildeS:badData', 'data.delta_range must be a cell array.');
  end
  if ~isfield(data,'Sdd_R') || ~iscell(data.Sdd_R) || ...
     ~isfield(data,'Scd')   || ~iscell(data.Scd)   || ...
     ~isfield(data,'Sdc')   || ~iscell(data.Sdc)
    error('solve_tildeS:badData', 'data.Sdd_R, data.Scd, data.Sdc must be cell arrays.');
  end

  sub  = data.sub;
  nSub = numel(sub);

  nC = data.primal.nC;
  nDeltaProd = data.nDeltaProd;

  % ADDED: validate scalar sizes
  if ~isnumeric(nC) || ~isscalar(nC) || ~isreal(nC) || ~isfinite(nC) || nC < 0 || fix(nC) ~= nC
    error('solve_tildeS:badData', 'data.primal.nC must be a nonnegative integer scalar.');
  end
  if ~isnumeric(nDeltaProd) || ~isscalar(nDeltaProd) || ~isreal(nDeltaProd) || ...
     ~isfinite(nDeltaProd) || nDeltaProd < 0 || fix(nDeltaProd) ~= nDeltaProd
    error('solve_tildeS:badData', 'data.nDeltaProd must be a nonnegative integer scalar.');
  end

  % ADDED: validate cell array lengths
  if numel(data.delta_range) ~= nSub || numel(data.Sdd_R) ~= nSub || ...
     numel(data.Scd) ~= nSub || numel(data.Sdc) ~= nSub
    error('solve_tildeS:badData', 'delta_range/Sdd_R/Scd/Sdc must have one entry per subdomain.');
  end

  % CHANGED: defaults are only via [] (not via missing args)
  if isempty(rc); rc = zeros(nC,1); end
  if isempty(rD); rD = zeros(nDeltaProd,1); end

  % ADDED: validate rc, rD
  if ~isnumeric(rc) || ~isreal(rc) || any(~isfinite(rc(:)))
    error('solve_tildeS:badInput', 'rc must be a real, finite numeric vector (or []).');
  end
  if ~isnumeric(rD) || ~isreal(rD) || any(~isfinite(rD(:)))
    error('solve_tildeS:badInput', 'rD must be a real, finite numeric vector (or []).');
  end

  rc = rc(:);
  rD = rD(:);

  if numel(rc) ~= nC
    error('solve_tildeS:badInput', 'rc must have length nC = %d (or be []).', nC);
  end
  if numel(rD) ~= nDeltaProd
    error('solve_tildeS:badInput', 'rD must have length nDeltaProd = %d (or be []).', nDeltaProd);
  end

  % ADDED: validate coarse factor if needed
  if nC > 0
    if ~isfield(data,'Kcc_R') || isempty(data.Kcc_R)
      error('solve_tildeS:badData', 'data.Kcc_R is required when nC > 0.');
    end
    Rcc = data.Kcc_R;
    if ~isnumeric(Rcc) || ~ismatrix(Rcc) || size(Rcc,1) ~= nC || size(Rcc,2) ~= nC
      error('solve_tildeS:badData', 'data.Kcc_R must be a %dx%d numeric matrix.', nC, nC);
    end
    if ~(istriu(Rcc) || istril(Rcc))
      error('solve_tildeS:badData', 'data.Kcc_R must be triangular (upper or lower).');
    end
  end

  % ADDED: validate indexing and operator dimensions; also ensure delta ranges are disjoint
  usedDelta = false(max(1,nDeltaProd), 1); % robust for nDeltaProd==0

  for i = 1:nSub
    % c_ids checks
    if ~isfield(sub(i),'c_ids')
      error('solve_tildeS:badData', 'data.sub(%d).c_ids is required.', i);
    end
    c_ids = sub(i).c_ids(:);

    if ~isempty(c_ids)
      if ~isnumeric(c_ids) || ~isreal(c_ids) || any(~isfinite(c_ids))
        error('solve_tildeS:badData', 'sub(%d).c_ids must be real, finite numeric indices.', i);
      end
      if any(fix(c_ids) ~= c_ids)
        error('solve_tildeS:badData', 'sub(%d).c_ids must be integer-valued.', i);
      end
      if nC == 0
        error('solve_tildeS:badData', 'sub(%d).c_ids must be empty when nC==0.', i);
      end
      if any(c_ids < 1) || any(c_ids > nC)
        error('solve_tildeS:badData', 'sub(%d).c_ids out of range 1..nC.', i);
      end
      if numel(unique(c_ids)) ~= numel(c_ids)
        error('solve_tildeS:badData', 'sub(%d).c_ids must be unique within the subdomain.', i);
      end
    end

    % delta_range checks
    rng = data.delta_range{i};
    if isempty(rng)
      localD = 0;
      rng = []; % normalize
    else
      if ~isnumeric(rng) || ~isreal(rng) || any(~isfinite(rng(:)))
        error('solve_tildeS:badData', 'delta_range{%d} must be real, finite numeric indices.', i);
      end
      rng = rng(:);
      if any(fix(rng) ~= rng)
        error('solve_tildeS:badData', 'delta_range{%d} must be integer-valued.', i);
      end
      if nDeltaProd == 0
        error('solve_tildeS:badData', 'delta_range{%d} must be empty when nDeltaProd==0.', i);
      end
      if any(rng < 1) || any(rng > nDeltaProd)
        error('solve_tildeS:badData', 'delta_range{%d} out of range 1..nDeltaProd.', i);
      end
      if any(usedDelta(rng))
        error('solve_tildeS:badData', 'delta_range{%d} overlaps another subdomain range.', i);
      end
      usedDelta(rng) = true;
      localD = numel(rng);
    end
    data.delta_range{i} = rng;

    % Operator validation
    Scd = data.Scd{i};
    Sdc = data.Sdc{i};

    if ~isnumeric(Scd) || ~ismatrix(Scd)
      error('solve_tildeS:badData', 'Scd{%d} must be a numeric matrix.', i);
    end
    if ~isnumeric(Sdc) || ~ismatrix(Sdc)
      error('solve_tildeS:badData', 'Sdc{%d} must be a numeric matrix.', i);
    end

    nCi = numel(c_ids);

    % FIXED: allow empty [] blocks for localD==0 (primal-only subdomains).
    if localD == 0
      % If provided non-empty, enforce correct "zero-column/zero-row" sizing.
      if ~isempty(Scd) && ~(size(Scd,1) == nCi && size(Scd,2) == 0)
        error('solve_tildeS:badData', ...
              'Scd{%d} must be empty or size (%d x 0) when delta_range is empty.', i, nCi);
      end
      if ~isempty(Sdc) && ~(size(Sdc,1) == 0 && size(Sdc,2) == nCi)
        error('solve_tildeS:badData', ...
              'Sdc{%d} must be empty or size (0 x %d) when delta_range is empty.', i, nCi);
      end
      % Sdd_R not required in this case
      continue;
    end

    % For localD > 0, enforce exact sizes
    if size(Scd,1) ~= nCi || size(Scd,2) ~= localD
      error('solve_tildeS:badData', ...
            'Scd{%d} must be size (%d x %d) to match c_ids and delta_range.', i, nCi, localD);
    end
    if size(Sdc,1) ~= localD || size(Sdc,2) ~= nCi
      error('solve_tildeS:badData', ...
            'Sdc{%d} must be size (%d x %d) to match delta_range and c_ids.', i, localD, nCi);
    end

    % Sdd_R checks
    R = data.Sdd_R{i};
    if ~isnumeric(R) || ~ismatrix(R) || size(R,1) ~= localD || size(R,2) ~= localD
      error('solve_tildeS:badData', 'Sdd_R{%d} must be size (%d x %d).', i, localD, localD);
    end
    if ~(istriu(R) || istril(R))
      error('solve_tildeS:badData', 'Sdd_R{%d} must be triangular (upper or lower).', i);
    end
  end

  % ------------------------------------------------------------
  % Step 1: yD = Sdd^{-1} rD (blockwise)
  % ------------------------------------------------------------
  yD = zeros(size(rD));
  yD_local = cell(nSub,1);

  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      yD_local{i} = zeros(0,1);
      continue;
    end
    ri = rD(rng);
    R = data.Sdd_R{i};
    yi = apply_fact_inv_(R, ri);
    yD(rng) = yi;
    yD_local{i} = yi;
  end

  % ------------------------------------------------------------
  % Step 2: coarse rhs
  % ------------------------------------------------------------
  rhs_c = rc;
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      continue;
    end
    yi = yD_local{i};
    if isempty(yi)
      continue;
    end
    rhs_c(c_ids) = rhs_c(c_ids) - data.Scd{i} * yi;
  end

  % ------------------------------------------------------------
  % Step 3: solve u_c = Kcc^{-1} rhs_c
  % ------------------------------------------------------------
  if nC > 0
    Rcc = data.Kcc_R;
    u_c = apply_fact_inv_(Rcc, rhs_c);
  else
    u_c = zeros(0,1);
  end

  % ------------------------------------------------------------
  % Step 4: local primal
  % ------------------------------------------------------------
  wc_local = cell(nSub,1);
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      wc_local{i} = zeros(0,1);
    else
      wc_local{i} = u_c(c_ids);
    end
  end

  % ------------------------------------------------------------
  % Step 5: wD correction
  % ------------------------------------------------------------
  wD = yD;
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    wci = wc_local{i};
    if isempty(wci)
      continue;
    end
    corr = data.Sdc{i} * wci;
    R = data.Sdd_R{i};
    z = apply_fact_inv_(R, corr);
    wD(rng) = wD(rng) - z;
  end

  out = struct();
  out.u_c = u_c;
  out.wc_local = wc_local;
  out.wD = wD;
end

% ============================================================
% Helper: apply inverse of SPD factor stored as triangular
% ============================================================
function x = apply_fact_inv_(T, b)
  if istriu(T) && ~istril(T)
    x = T \ (T' \ b);
  elseif istril(T) && ~istriu(T)
    x = T' \ (T \ b);
  else
    x = T \ (T' \ b);
  end
end