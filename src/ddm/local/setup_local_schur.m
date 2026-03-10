function sub = setup_local_schur(sub, opts)
%SETUP_LOCAL_SCHUR Prepare local Schur-complement data on each subdomain.
% Thesis link: Chapter 4.2.2–4.2.4 (interior elimination and interface reduction).
% The routine factorizes `K_II`, builds the reduced interface RHS, and
% optionally assembles the explicit local Schur complement.
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

  % keep signature strict but clearer
  if nargin < 1 || nargin > 2
    error('setup_local_schur: expected inputs (sub[,opts]).');
  end

  % sub must be struct
  if ~isstruct(sub)
    error('setup_local_schur: sub must be a struct array.');
  end

  % robust opts handling
  if nargin == 1
    opts = struct();
  end
  % opts must be struct (avoid generic isfield crash)
  if ~isstruct(opts)
    error('setup_local_schur: opts must be a struct.');
  end
  if ~isfield(opts,'assemble_S')
    opts.assemble_S = false;
  end
  % assemble_S must be scalar flag
  if ischar(opts.assemble_S) || ~isscalar(opts.assemble_S)
    error('setup_local_schur: opts.assemble_S must be a scalar logical/numeric flag.');
  end
  opts.assemble_S = logical(opts.assemble_S);

  %if not assembling S, remove any pre-existing S safely for whole array
  if ~opts.assemble_S && isfield(sub,'S')
    sub = rmfield(sub,'S');
  end

  nSub = numel(sub);

  % required field list (block matrices + RHS)
  req = {'K_II','K_Ig','K_gI','K_gg','f_I','f_g'};

  for i = 1:nSub
    % previous version checked only K_II; now require all needed fields
    for k = 1:numel(req)
      if ~isfield(sub(i), req{k})
        error('setup_local_schur: sub(%d) missing field %s (run extract_subdomain_blocks first).', i, req{k});
      end
    end

    % ---- fetch local copies (helps validate sizes cleanly)
    K_II = sub(i).K_II;
    K_Ig = sub(i).K_Ig;
    K_gI = sub(i).K_gI;
    K_gg = sub(i).K_gg;
    f_I  = sub(i).f_I;
    f_g  = sub(i).f_g;

    %reject invalid/ambiguous types; require numeric real finite
    validate_numeric_real_finite_('K_II', K_II, i);  % rejects char/logical/complex/NaN/Inf
    validate_numeric_real_finite_('K_Ig', K_Ig, i);
    validate_numeric_real_finite_('K_gI', K_gI, i);
    validate_numeric_real_finite_('K_gg', K_gg, i);
    validate_numeric_real_finite_('f_I',  f_I,  i);
    validate_numeric_real_finite_('f_g',  f_g,  i);

    % block size checks
    nI = size(K_II,1);
    if size(K_II,2) ~= nI
      error('setup_local_schur: sub(%d).K_II must be square.', i);
    end

    nG = size(K_gg,1);
    if size(K_gg,2) ~= nG
      error('setup_local_schur: sub(%d).K_gg must be square.', i);
    end

    if ~isequal(size(K_Ig), [nI, nG])
      error('setup_local_schur: sub(%d).K_Ig must be size %dx%d.', i, nI, nG);
    end
    if ~isequal(size(K_gI), [nG, nI])
      error('setup_local_schur: sub(%d).K_gI must be size %dx%d.', i, nG, nI);
    end

    %RHS must be vectors of correct lengths; normalize to column
    if ~(isvector(f_I) || isempty(f_I))
      error('setup_local_schur: sub(%d).f_I must be a vector.', i);
    end
    if ~(isvector(f_g) || isempty(f_g))
      error('setup_local_schur: sub(%d).f_g must be a vector.', i);
    end
    f_I = f_I(:);
    f_g = f_g(:);
    if numel(f_I) ~= nI
      error('setup_local_schur: sub(%d).f_I must have length %d.', i, nI);
    end
    if numel(f_g) ~= nG
      error('setup_local_schur: sub(%d).f_g must have length %d.', i, nG);
    end

    %nI==0 branch uses validated/column f_g
    if nI == 0
      sub(i).R_II = [];
      sub(i).g = f_g;

      if opts.assemble_S
        sub(i).S = K_gg;
      end
      continue;
    end

    % Cholesky factorization (SPD).
    % chol returns upper triangular R such that R'*R = K_II.
    [R, pflag] = chol(K_II);
    if pflag ~= 0
      error('setup_local_schur: chol failed for sub(%d).K_II (matrix not SPD?).', i);
    end
    sub(i).R_II = R;

    % Reduced RHS: g = f_g - K_gI * (K_II^{-1} f_I)
    z = R \ (R' \ f_I); % use validated/column f_I
    sub(i).g = f_g - K_gI * z; %use validated/column f_g

    if opts.assemble_S
      X = R \ (R' \ K_Ig);        % use validated K_Ig
      sub(i).S = K_gg - K_gI * X; % use validated K_gg/K_gI
    end
  end
end

% =========================
% Helpers (file-private)
% =========================

function validate_numeric_real_finite_(name, v, i)
% shared validation (Rule 7.1 hardening target)
  if ischar(v) || islogical(v) || ~isnumeric(v) || ~isreal(v)
    error('setup_local_schur: sub(%d).%s must be numeric real (reject char/logical).', i, name);
  end

  if issparse(v)
    nz = nonzeros(v);
    if any(~isfinite(nz))
      error('setup_local_schur: sub(%d).%s must be finite (no NaN/Inf).', i, name);
    end
  else
    if any(~isfinite(v(:)))
      error('setup_local_schur: sub(%d).%s must be finite (no NaN/Inf).', i, name);
    end
  end
end