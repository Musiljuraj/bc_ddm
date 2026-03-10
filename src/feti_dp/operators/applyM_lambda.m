% ============================================================
% File: src/feti_dp/operators/applyM_lambda.m
% ============================================================
function z = applyM_lambda(r, data)
%APPLYM_LAMBDA Apply the baseline FETI-DP preconditioner in multiplier space.
% Thesis link: Chapter 5.3 and Chapter 6.2 (preconditioned FETI-DP iteration).
% The routine applies scaling, local `S_{\Delta\Delta}` blocks, and `B_d`
% in the matrix-free preconditioner action `M_lambda^{-1} r`.

  % ----------------------------
  % Input validation / hardening
  % ----------------------------

  %  : strict arg-count check
  if nargin ~= 2
    error('applyM_lambda:badNargin', 'applyM_lambda expects exactly 2 inputs: (r, data).');
  end

  %  : validate r type/value
  if ~isnumeric(r) || islogical(r) || ~isreal(r)
    error('applyM_lambda:badRType', 'r must be a real, non-logical numeric vector.');
  end
  if isempty(r) || ~isvector(r)
    error('applyM_lambda:badRShape', 'r must be a non-empty vector.');
  end
  if any(~isfinite(r(:)))
    error('applyM_lambda:badRFinite', 'r must contain only finite values (no NaN/Inf).');
  end

  %  : validate data struct + required fields
  if ~isstruct(data)
    error('applyM_lambda:badDataType', 'data must be a struct.');
  end
  req = {'BdT','Bd','DeltaWeights','sub','delta_range','Sdd'};
  for k = 1:numel(req)
    if ~isfield(data, req{k})
      error('applyM_lambda:missingField', 'data.%s is required.', req{k});
    end
  end

  %  : validate operator containers
  if ~isnumeric(data.BdT) || ~isreal(data.BdT)
    error('applyM_lambda:badBdT', 'data.BdT must be a real numeric matrix.');
  end
  if ~isnumeric(data.Bd)  || ~isreal(data.Bd)
    error('applyM_lambda:badBd', 'data.Bd must be a real numeric matrix.');
  end
  if ~isnumeric(data.DeltaWeights) || islogical(data.DeltaWeights) || ~isreal(data.DeltaWeights)
    error('applyM_lambda:badDeltaWeightsType', 'data.DeltaWeights must be a real numeric vector.');
  end
  if ~isvector(data.DeltaWeights) || isempty(data.DeltaWeights)
    error('applyM_lambda:badDeltaWeightsShape', 'data.DeltaWeights must be a non-empty vector.');
  end
  if any(~isfinite(data.DeltaWeights(:)))
    error('applyM_lambda:badDeltaWeightsFinite', 'data.DeltaWeights must be finite (no NaN/Inf).');
  end
  if ~iscell(data.delta_range) || ~iscell(data.Sdd)
    error('applyM_lambda:badCells', 'data.delta_range and data.Sdd must be cell arrays.');
  end

  %  : normalize r to column *after* validation
  r = r(:);

  nSub = numel(data.sub);

  %  : length consistency between sub, delta_range, Sdd
  if numel(data.delta_range) ~= nSub || numel(data.Sdd) ~= nSub
    error('applyM_lambda:badSubCount', ...
          'data.delta_range and data.Sdd must have the same length as data.sub.');
  end

  %  : dimension consistency checks for BdT and r
  if ndims(data.BdT) ~= 2 || size(data.BdT,2) ~= numel(r)
    error('applyM_lambda:dimMismatchBdT', ...
          'Size mismatch: size(BdT,2) must equal numel(r).');
  end

  % 1) t = Bd^T r
  t = data.BdT * r;                          % packed Delta

  %  : DeltaWeights must match packed-Delta size
  if numel(data.DeltaWeights) ~= numel(t)
    error('applyM_lambda:dimMismatchDeltaWeights', ...
          'Size mismatch: numel(DeltaWeights) must equal size(BdT,1).');
  end

  %  : Bd must map packed-Delta -> lambda space
  if ndims(data.Bd) ~= 2 || size(data.Bd,2) ~= numel(t)
    error('applyM_lambda:dimMismatchBd', ...
          'Size mismatch: size(Bd,2) must equal size(BdT,1).');
  end

  % 2) t = D t
  t = data.DeltaWeights(:) .* t;

  % 3) u = Sdd * t (blockwise multiply, no solves)
  u = zeros(size(t));

  %  : range checks (integer, in-bounds, no duplicates/overlaps)
  used = false(numel(t), 1);
  for i = 1:nSub
    rng = data.delta_range{i};

    if isempty(rng)
      continue;
    end

    if ~isnumeric(rng) || islogical(rng) || ~isreal(rng) || ~isvector(rng)
      error('applyM_lambda:badRangeType', 'delta_range{%d} must be a real numeric index vector.', i);
    end
    rng = rng(:);

    if any(~isfinite(rng))
      error('applyM_lambda:badRangeFinite', 'delta_range{%d} contains NaN/Inf.', i);
    end
    if any(abs(rng - round(rng)) > 0)
      error('applyM_lambda:badRangeInteger', 'delta_range{%d} must contain integer-valued indices.', i);
    end
    rng = round(rng);

    if any(rng < 1) || any(rng > numel(t))
      error('applyM_lambda:badRangeBounds', ...
            'delta_range{%d} contains out-of-range indices (valid: 1..%d).', i, numel(t));
    end
    if numel(unique(rng)) ~= numel(rng)
      error('applyM_lambda:badRangeDuplicates', 'delta_range{%d} contains duplicate indices.', i);
    end
    if any(used(rng))
      error('applyM_lambda:badRangeOverlap', ...
            'delta_range{%d} overlaps with another range (order-dependent overwrite).', i);
    end
    used(rng) = true;

    %  : validate local block operator Sdd{i}
    A = data.Sdd{i};
    if ~isnumeric(A) || islogical(A) || ~isreal(A) || ndims(A) ~= 2
      error('applyM_lambda:badSddType', 'Sdd{%d} must be a real numeric matrix.', i);
    end
    if any(~isfinite(A(:)))
      error('applyM_lambda:badSddFinite', 'Sdd{%d} must contain only finite values (no NaN/Inf).', i);
    end
    m = numel(rng);
    if size(A,1) ~= m || size(A,2) ~= m
      error('applyM_lambda:badSddSize', ...
            'Sdd{%d} must be %d-by-%d to match numel(delta_range{%d}).', i, m, m, i);
    end

    ti = t(rng);
    u(rng) = A * ti;
  end

  % 4) u = D u
  u = data.DeltaWeights(:) .* u;

  % 5) z = Bd u
  z = data.Bd * u;
end