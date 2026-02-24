function [B, meta] = build_jump_operator_B(prod)
%BUILD_JUMP_OPERATOR_B  Construct the jump/constraint operator B enforcing continuity of duplicated interface DOFs.
%
% Link to thesis:
%   Chapter 3.3.3 (jump operator), equations (3.42)–(3.44).
%
% In the product interface space W, each physical interface DOF (assembled index "hat")
% may have k>=2 local copies (product indices). Continuity requires that all copies match.
% We enforce this by constraints:
%   w(p1) - w(pj) = 0,  j=2..k,
% where p1 is a chosen reference copy.
%
% Inputs:
%   prod : product-space bookkeeping from build_product_interface (needs .hat2prod and .nProd/.nHat).
%
% Outputs:
%   B    : sparse constraint matrix of size (m x nProd), with entries in {-1,0,+1}.
%   meta : (optional) struct with fields:
%          .row_hat : (m x 1) assembled interface index associated with each constraint row

  if nargin ~= 1
    error('build_jump_operator_B:nargin', 'Expected exactly 1 input (prod).'); % CHANGED
  end
  if ~isstruct(prod) || ~isfield(prod,'hat2prod') || ~isfield(prod,'nProd') || ~isfield(prod,'nHat')
    error('build_jump_operator_B:badProd', ...
          'prod must contain hat2prod, nProd, nHat (run build_product_interface).'); % CHANGED
  end

  nHat  = prod.nHat;
  nProd = prod.nProd;

  % Basic scalar validations (solid standard).                                  % ADDED
  if ~isnumeric(nHat) || ~isscalar(nHat) || ~isreal(nHat) || ~isfinite(nHat) || nHat ~= round(nHat) || nHat < 0
    error('build_jump_operator_B:bad_nHat', 'nHat must be a real, finite, nonnegative integer scalar.'); % ADDED
  end
  if ~isnumeric(nProd) || ~isscalar(nProd) || ~isreal(nProd) || ~isfinite(nProd) || nProd ~= round(nProd) || nProd < 0
    error('build_jump_operator_B:bad_nProd', 'nProd must be a real, finite, nonnegative integer scalar.'); % ADDED
  end
  if ~iscell(prod.hat2prod)
    error('build_jump_operator_B:bad_hat2prod', 'hat2prod must be a cell array.'); % ADDED
  end
  if numel(prod.hat2prod) < nHat
    error('build_jump_operator_B:bad_hat2prod_len', 'hat2prod must have at least nHat entries.'); % ADDED
  end

  % Count constraints: sum over hat DOFs of (k-1).
  m = 0;
  for h = 1:nHat
    copies = validate_copies_(prod.hat2prod{h}, h, nProd); % ADDED
    k = numel(copies);
    if k > 1
      m = m + (k-1);
    end
  end

  % Triplets for sparse B: each row has exactly two nonzeros (+1 and -1).
  I = zeros(2*m, 1);
  J = zeros(2*m, 1);
  V = zeros(2*m, 1);

  row_hat = zeros(m, 1);

  row = 1;
  ptr = 1;

  % Cycle all w-hat elements, take all copies and create constraints for them.
  for h = 1:nHat
    copies = validate_copies_(prod.hat2prod{h}, h, nProd); % CHANGED (was raw access)
    k = numel(copies);
    if k <= 1
      continue;
    end

    p1 = copies(1);
    for j = 2:k
      pj = copies(j);

      % Constraint row: w(p1) - w(pj) = 0
      I(ptr)   = row; J(ptr)   = p1; V(ptr)   =  1; ptr = ptr + 1;
      I(ptr)   = row; J(ptr)   = pj; V(ptr)   = -1; ptr = ptr + 1;

      row_hat(row) = h;
      row = row + 1;
    end
  end

  % Internal consistency (cheap, avoids silent size mismatch).                   % ADDED
  if row ~= (m + 1) || ptr ~= (2*m + 1)
    error('build_jump_operator_B:internal', 'Internal assembly mismatch (row/ptr counters).'); % ADDED
  end

  B = sparse(I, J, V, m, nProd);

  if nargout >= 2
    meta = struct();
    meta.row_hat = row_hat;
  end
end

% -------------------------------------------------------------------------
function copies = validate_copies_(copies, h, nProd)
%validate_copies_  Minimal validation for a hat->prod index list.            % ADDED
%
% Requirements (solid standard):
%   - nonempty numeric vector
%   - real, finite
%   - integer-valued indices in [1..nProd]
%   - no duplicates (avoids zero constraint rows)

  if isempty(copies)
    error('build_jump_operator_B:emptyCopies', 'hat2prod{%d} is empty.', h); % ADDED
  end
  if ~isnumeric(copies)
    error('build_jump_operator_B:badCopiesType', 'hat2prod{%d} must be numeric indices.', h); % ADDED
  end

  copies = copies(:);

  if ~isreal(copies) || any(~isfinite(copies))
    error('build_jump_operator_B:badCopiesValues', 'hat2prod{%d} must be real and finite.', h); % ADDED
  end
  if any(copies ~= round(copies))
    error('build_jump_operator_B:badCopiesInt', 'hat2prod{%d} must be integer-valued.', h); % ADDED
  end
  if any(copies < 1) || any(copies > nProd)
    error('build_jump_operator_B:badCopiesRange', 'hat2prod{%d} contains indices out of range 1..nProd.', h); % ADDED
  end
  if numel(unique(copies)) ~= numel(copies)
    error('build_jump_operator_B:dupCopies', 'hat2prod{%d} contains duplicate indices.', h); % ADDED
  end
end