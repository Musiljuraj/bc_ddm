% ============================================================
% File: src/feti_dp/solve/reconstruct_fetidp_solution.m
% ============================================================
function [w_c, w_d, u_free, diag] = reconstruct_fetidp_solution(lambda, data)
%RECONSTRUCT_FETIDP_SOLUTION  Recover interface traces and full free-DOF solution from lambda.
%
% Outputs:
%   w_c   : coarse vector u_c (assembled primal dofs, size nC)
%   w_d   : packed Delta vector (size nDeltaProd)
%   u_free: global solution in FREE-DOF numbering (size nFree x 1)
%   diag  : diagnostics (constraint residual norm, etc.)

  narginchk(2,2); % ADDED

  % ---- basic input validation (light hardening) ----
  if ~(isnumeric(lambda) && isreal(lambda)) % ADDED
    error('reconstruct_fetidp_solution:lambdaType', 'lambda must be a real numeric vector'); % ADDED
  end % ADDED
  lambda = lambda(:); % CHANGED (moved after type check)
  if any(~isfinite(lambda)) % ADDED
    error('reconstruct_fetidp_solution:lambdaFinite', 'lambda must contain only finite values'); % ADDED
  end % ADDED

  if ~isstruct(data) % ADDED
    error('reconstruct_fetidp_solution:dataType', 'data must be a struct'); % ADDED
  end % ADDED

  req_fields = {'sub','nDeltaProd','delta_range','BdT','Sdd_R','primal','hc','Scd','Sdc','free','Bd'}; % ADDED
  for k = 1:numel(req_fields) % ADDED
    if ~isfield(data, req_fields{k}) % ADDED
      error('reconstruct_fetidp_solution:missingField', 'data is missing required field: %s', req_fields{k}); % ADDED
    end % ADDED
  end % ADDED
  if ~isfield(data.primal,'nC') % ADDED
    error('reconstruct_fetidp_solution:missingField', 'data.primal is missing required field: nC'); % ADDED
  end % ADDED

  if ~(isscalar(data.nDeltaProd) && isnumeric(data.nDeltaProd) && isreal(data.nDeltaProd) && isfinite(data.nDeltaProd)) % ADDED
    error('reconstruct_fetidp_solution:nDeltaProd', 'data.nDeltaProd must be a real finite scalar'); % ADDED
  end % ADDED
  if data.nDeltaProd < 0 || fix(data.nDeltaProd) ~= data.nDeltaProd % ADDED
    error('reconstruct_fetidp_solution:nDeltaProdInt', 'data.nDeltaProd must be a nonnegative integer'); % ADDED
  end % ADDED

  if ~iscell(data.delta_range) % ADDED
    error('reconstruct_fetidp_solution:deltaRangeType', 'data.delta_range must be a cell array'); % ADDED
  end % ADDED
  if ~iscell(data.Sdd_R) || ~iscell(data.Scd) || ~iscell(data.Sdc) % ADDED
    error('reconstruct_fetidp_solution:cellFields', 'data.Sdd_R, data.Scd, data.Sdc must be cell arrays'); % ADDED
  end % ADDED

  if ~(isnumeric(data.BdT) && isreal(data.BdT)) % ADDED
    error('reconstruct_fetidp_solution:BdTType', 'data.BdT must be a real numeric matrix'); % ADDED
  end % ADDED

  if ~(isnumeric(data.Bd) && isreal(data.Bd)) % ADDED
    error('reconstruct_fetidp_solution:BdType', 'data.Bd must be a real numeric matrix'); % ADDED
  end % ADDED
  if size(data.Bd,2) ~= data.nDeltaProd % ADDED
    error('reconstruct_fetidp_solution:BdSize', 'size(data.Bd,2) must equal data.nDeltaProd'); % ADDED
  end % ADDED

  if ~(isvector(data.free) && isnumeric(data.free)) % ADDED
    error('reconstruct_fetidp_solution:free', 'data.free must be a numeric vector'); % ADDED
  end % ADDED
  nFree = numel(data.free);

  sub = data.sub;
  if ~isstruct(sub) % ADDED
    error('reconstruct_fetidp_solution:subType', 'data.sub must be a struct array'); % ADDED
  end % ADDED
  nSub = numel(sub);

  if numel(data.delta_range) ~= nSub % ADDED
    error('reconstruct_fetidp_solution:deltaRangeLen', 'numel(data.delta_range) must equal numel(data.sub)'); % ADDED
  end % ADDED
  if numel(data.Sdd_R) ~= nSub || numel(data.Scd) ~= nSub || numel(data.Sdc) ~= nSub % ADDED
    error('reconstruct_fetidp_solution:cellLen', 'data.Sdd_R / Scd / Sdc must have one entry per subdomain'); % ADDED
  end % ADDED

  if size(data.BdT,1) ~= data.nDeltaProd % ADDED
    error('reconstruct_fetidp_solution:BdTSize1', 'size(data.BdT,1) must equal data.nDeltaProd'); % ADDED
  end % ADDED
  if size(data.BdT,2) ~= numel(lambda) % ADDED
    error('reconstruct_fetidp_solution:BdTSize2', 'size(data.BdT,2) must equal numel(lambda)'); % ADDED
  end % ADDED

  % Validate per-subdomain index consistency that is relied upon later
  nC = data.primal.nC; % ADDED
  for i = 1:nSub % ADDED
    sub_req = {'g','idx_c','idx_d','c_ids','glob_G'}; % ADDED
    for k = 1:numel(sub_req) % ADDED
      if ~isfield(sub(i), sub_req{k}) % ADDED
        error('reconstruct_fetidp_solution:missingSubField', 'sub(%d) is missing required field: %s', i, sub_req{k}); % ADDED
      end % ADDED
    end % ADDED

    rng = data.delta_range{i}; % ADDED
    if ~isempty(rng) % ADDED
      validate_index_vec_(rng, data.nDeltaProd, sprintf('data.delta_range{%d}', i)); % ADDED
    end % ADDED

    idx_c = sub(i).idx_c(:); % ADDED
    idx_d = sub(i).idx_d(:); % ADDED
    if ~isempty(idx_c), validate_index_vec_(idx_c, inf, sprintf('sub(%d).idx_c', i)); end % ADDED
    if ~isempty(idx_d), validate_index_vec_(idx_d, inf, sprintf('sub(%d).idx_d', i)); end % ADDED

    c_ids = sub(i).c_ids(:); % ADDED
    if ~isempty(c_ids) % ADDED
      validate_index_vec_(c_ids, nC, sprintf('sub(%d).c_ids', i)); % ADDED
    end % ADDED

    globG = sub(i).glob_G(:); % ADDED
    validate_index_vec_(globG, nFree, sprintf('sub(%d).glob_G', i)); % ADDED

    nG = numel(globG); % ADDED
    if ~isempty([idx_c; idx_d]) % ADDED
      mx = max([idx_c; idx_d]); % ADDED
      if mx > nG % ADDED
        error('reconstruct_fetidp_solution:gammaIndexRange', ...
              'sub(%d) idx_c/idx_d exceed gamma length implied by numel(glob_G)', i); % ADDED
      end % ADDED
    end % ADDED

    if ~isempty(rng) && ~isempty(idx_d) % ADDED
      if numel(rng) ~= numel(idx_d) % ADDED
        error('reconstruct_fetidp_solution:deltaMapLen', ...
              'sub(%d) mismatch: numel(delta_range) ~= numel(idx_d)', i); % ADDED
      end % ADDED
    end % ADDED
    if isempty(rng) && ~isempty(idx_d) % ADDED
      error('reconstruct_fetidp_solution:deltaMapEmpty', ...
            'sub(%d) has nonempty idx_d but empty delta_range', i); % ADDED
    end % ADDED

    if ~isempty(idx_c) && ~isempty(c_ids) % ADDED
      if numel(idx_c) ~= numel(c_ids) % ADDED
        error('reconstruct_fetidp_solution:coarseMapLen', ...
              'sub(%d) mismatch: numel(idx_c) ~= numel(c_ids)', i); % ADDED
      end % ADDED
    end % ADDED

    % Check local factor sizes (Sdd_R) against delta block length
    if ~isempty(rng) % ADDED
      R = data.Sdd_R{i}; % ADDED
      if ~(isnumeric(R) && isreal(R)) % ADDED
        error('reconstruct_fetidp_solution:SddRType', 'data.Sdd_R{%d} must be real numeric', i); % ADDED
      end % ADDED
      if size(R,1) ~= size(R,2) || size(R,1) ~= numel(rng) % ADDED
        error('reconstruct_fetidp_solution:SddRSize', ...
              'data.Sdd_R{%d} must be square of size numel(delta_range{%d})', i, i); % ADDED
      end % ADDED

      % Light checks for coupling operator sizes where used
      if ~isempty(c_ids) % ADDED
        Scd = data.Scd{i}; % ADDED
        Sdc = data.Sdc{i}; % ADDED
        if ~(isnumeric(Scd) && isreal(Scd)) || ~(isnumeric(Sdc) && isreal(Sdc)) % ADDED
          error('reconstruct_fetidp_solution:couplingType', 'data.Scd{%d}/Sdc{%d} must be real numeric', i, i); % ADDED
        end % ADDED
        if size(Scd,2) ~= numel(rng) % ADDED
          error('reconstruct_fetidp_solution:ScdSize', 'size(Scd{%d},2) must equal numel(delta_range{%d})', i, i); % ADDED
        end % ADDED
        if size(Sdc,1) ~= numel(rng) || size(Sdc,2) ~= numel(c_ids) % ADDED
          error('reconstruct_fetidp_solution:SdcSize', ...
                'size(Sdc{%d}) must be [numel(delta_range{%d}) x numel(c_ids)]', i, i); % ADDED
        end % ADDED
      end % ADDED
    end % ADDED

    % Check g length supports idx_d extraction (and gamma length implied by glob_G)
    gi = sub(i).g(:); % ADDED
    if ~(isnumeric(gi) && isreal(gi) && all(isfinite(gi))) % ADDED
      error('reconstruct_fetidp_solution:gType', 'sub(%d).g must be real finite numeric', i); % ADDED
    end % ADDED
    if numel(gi) < nG % ADDED
      error('reconstruct_fetidp_solution:gLen', 'sub(%d).g length is smaller than numel(sub(%d).glob_G)', i, i); % ADDED
    end % ADDED

    % Validate optional interior scatter indices
    if isfield(sub(i), 'glob_I') && ~isempty(sub(i).glob_I) % ADDED
      validate_index_vec_(sub(i).glob_I(:), nFree, sprintf('sub(%d).glob_I', i)); % ADDED
    end % ADDED
  end % ADDED

  if ~(isnumeric(data.hc) && isreal(data.hc)) % ADDED
    error('reconstruct_fetidp_solution:hcType', 'data.hc must be a real numeric vector'); % ADDED
  end % ADDED
  if numel(data.hc) ~= nC % ADDED
    error('reconstruct_fetidp_solution:hcSize', 'numel(data.hc) must equal data.primal.nC'); % ADDED
  end % ADDED
  if nC > 0 % ADDED
    if ~isfield(data,'Kcc_R') % ADDED
      error('reconstruct_fetidp_solution:missingField', 'data is missing required field: Kcc_R (needed when nC>0)'); % ADDED
    end % ADDED
    Rcc = data.Kcc_R; % ADDED
    if ~(isnumeric(Rcc) && isreal(Rcc)) % ADDED
      error('reconstruct_fetidp_solution:KccRType', 'data.Kcc_R must be real numeric'); % ADDED
    end % ADDED
    if size(Rcc,1) ~= size(Rcc,2) || size(Rcc,1) ~= nC % ADDED
      error('reconstruct_fetidp_solution:KccRSize', 'data.Kcc_R must be square of size nC'); % ADDED
    end % ADDED
  end % ADDED

  % ---- original computation ----

  % Build gD packed and gC (data.hC) assembled into coarse ids using hc already computed.
  gD = zeros(data.nDeltaProd,1);
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    gi = sub(i).g(:);
    gD(rng) = gi(sub(i).idx_d(:));
  end

  % C*lambda contribution to coarse rhs:
  % Compute y = Sdd^{-1} (Bd^T lambda) and assemble -Scd*y into coarse rhs
  t = data.BdT * lambda;
  y_local = cell(nSub,1);
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      y_local{i} = zeros(0,1);
      continue;
    end
    ti = t(rng);
    R = data.Sdd_R{i};
    yi = R \ (R' \ ti);
    y_local{i} = yi;
  end

  rhs_c = data.hc;
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      continue;
    end
    yi = y_local{i};
    if isempty(yi)
      continue;
    end
    rhs_c(c_ids) = rhs_c(c_ids) + data.Scd{i} * yi; % sign consistent with -Bd^T lambda inside delta equation
  end

  % Solve u_c
  if data.primal.nC > 0
    Rcc = data.Kcc_R;
    u_c = Rcc \ (Rcc' \ rhs_c);
  else
    u_c = zeros(0,1);
  end

  % Local primal traces wci
  wc_local = cell(nSub,1);
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      wc_local{i} = zeros(0,1);
    else
      wc_local{i} = u_c(c_ids);
    end
  end

  % Compute delta trace:
  % wD = Sdd^{-1}( gD - Sdc*w_c - Bd^T lambda )
  rhsD = gD - t;
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    wci = wc_local{i};
    if ~isempty(wci)
      rhsD(rng) = rhsD(rng) - data.Sdc{i} * wci;
    end
  end

  wD = zeros(size(rhsD));
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    R = data.Sdd_R{i};
    wD(rng) = R \ (R' \ rhsD(rng));
  end

  % Reconstruct local gamma vectors and interior u_I via back-substitution
  u_free_accum = zeros(nFree,1);
  u_free_count = zeros(nFree,1);

  for i = 1:nSub
    % local gamma in original sub(i).gamma ordering
    globG = sub(i).glob_G(:);                 % ADDED (reuse)
    nG = numel(globG);                        % FIXED (was numel(idx_c)+numel(idx_d))
    w_gamma = zeros(nG,1);

    % fill c and d parts
    w_gamma(sub(i).idx_c(:)) = wc_local{i};
    if ~isempty(sub(i).idx_d)
      w_gamma(sub(i).idx_d(:)) = wD(data.delta_range{i});
    end

    % Recover interior for this subdomain (if any)
    if isfield(sub(i), 'R_II') && ~isempty(sub(i).R_II)
      % light consistency checks for interior solve inputs % ADDED
      if ~isfield(sub(i),'f_I') || ~isfield(sub(i),'K_Ig') % ADDED
        error('reconstruct_fetidp_solution:missingInteriorFields', ...
              'sub(%d) has R_II but is missing f_I and/or K_Ig', i); % ADDED
      end % ADDED
      RII = sub(i).R_II; % ADDED
      if ~(isnumeric(RII) && isreal(RII)) % ADDED
        error('reconstruct_fetidp_solution:RIIType', 'sub(%d).R_II must be real numeric', i); % ADDED
      end % ADDED
      if size(RII,1) ~= size(RII,2) % ADDED
        error('reconstruct_fetidp_solution:RIISquare', 'sub(%d).R_II must be square', i); % ADDED
      end % ADDED
      fI = sub(i).f_I(:); % ADDED
      KIg = sub(i).K_Ig;  % ADDED
      if size(KIg,2) ~= nG % ADDED
        error('reconstruct_fetidp_solution:KIgSize', 'sub(%d).K_Ig must have %d columns (nG)', i, nG); % ADDED
      end % ADDED
      if size(KIg,1) ~= numel(fI) || size(RII,1) ~= numel(fI) % ADDED
        error('reconstruct_fetidp_solution:InteriorSize', 'sub(%d) interior sizes inconsistent (f_I, K_Ig, R_II)', i); % ADDED
      end % ADDED

      tmp = fI - KIg * w_gamma; % CHANGED (use locals)
      uI = RII \ (RII' \ tmp);  % CHANGED (use locals)
    else
      uI = zeros(0,1);
    end

    % Scatter to global FREE-DOF vector by averaging duplicates
    u_free_accum(globG) = u_free_accum(globG) + w_gamma;
    u_free_count(globG) = u_free_count(globG) + 1;

    % Interior free dofs:
    if isfield(sub(i), 'glob_I') && ~isempty(sub(i).glob_I)
      globI = sub(i).glob_I(:);
      u_free_accum(globI) = u_free_accum(globI) + uI;
      u_free_count(globI) = u_free_count(globI) + 1;
    end
  end

  % Final averaging for duplicates
  mask = (u_free_count > 0);
  u_free = zeros(nFree,1);
  u_free(mask) = u_free_accum(mask) ./ u_free_count(mask);

  % Outputs
  w_c = u_c;
  w_d = wD;

  diag = struct();
  diag.constraint_norm = norm(data.Bd * wD, 2); % CHANGED (avoid full())
end

% ---- local helpers ----
function validate_index_vec_(idx, nMax, name)
  % ADDED
  if ~(isnumeric(idx) && isreal(idx))
    error('reconstruct_fetidp_solution:indexType', '%s must be a real numeric index vector', name);
  end
  idx = idx(:);
  if isempty(idx)
    return;
  end
  if any(~isfinite(idx))
    error('reconstruct_fetidp_solution:indexFinite', '%s must be finite', name);
  end
  if any(idx < 1) || any(fix(idx) ~= idx)
    error('reconstruct_fetidp_solution:indexInt', '%s must contain positive integer indices', name);
  end
  if isfinite(nMax) && any(idx > nMax)
    error('reconstruct_fetidp_solution:indexRange', '%s indices out of range (max allowed %d)', name, nMax);
  end
end