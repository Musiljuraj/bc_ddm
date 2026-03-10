% ============================================================
% File: src/common/spectra/assemble_from_apply.m
% ============================================================
function Mat = assemble_from_apply(applyFun, n, opts)
%ASSEMBLE_FROM_APPLY Assemble an explicit matrix from a matrix-free callback.
% Thesis link: Chapter 6.3 (spectral analysis in the matrix-free setting).
% For small problems, the routine builds the matrix column-by-column from
% repeated operator applications to standard basis vectors.
%
% Inputs:
%   applyFun : function handle, y = applyFun(x), where x,y are length-n vectors.
%              (In this thesis, applyFun will wrap applyA_* or applyM_* routines.)
%   n        : positive integer, dimension of the operator domain/range.
%   opts     : optional struct with fields:
%              - verbose         : if true, print coarse progress information.
%              - progress_every  : integer >= 1, print every k columns if verbose.
%              - force_full      : if true (default), convert outputs to full.
%
% Output:
%   Mat : explicit n-by-n dense matrix whose j-th column is applyFun(e_j).
%
  % ----------------------------
  % Input validation
  % ----------------------------
  if nargin < 2 || nargin > 3
    error('assemble_from_apply:InvalidNargin: expected inputs (applyFun, n [, opts]).');
  end

  if ~isa(applyFun, 'function_handle')
    error('assemble_from_apply:InvalidApplyFun: applyFun must be a function handle.');
  end

  if ~(isscalar(n) && isnumeric(n) && isreal(n) && isfinite(n) && n == floor(n) && n >= 1)
    error('assemble_from_apply:InvalidN: n must be a positive integer scalar.');
  end

  if nargin < 3 || isempty(opts)
    opts = struct();
  end
  if ~isstruct(opts)
    error('assemble_from_apply:InvalidOpts: opts must be a struct (or []).');
  end

  % ----------------------------
  % Options (defaults)
  % ----------------------------
  if ~isfield(opts, 'verbose') || isempty(opts.verbose)
    opts.verbose = false;
  end
  if ~islogical(opts.verbose) || ~isscalar(opts.verbose)
    error('assemble_from_apply:InvalidVerbose: opts.verbose must be a logical scalar.');
  end

  if ~isfield(opts, 'progress_every') || isempty(opts.progress_every)
    % Default: print about ~10 updates if n is moderate, otherwise at least 1.
    opts.progress_every = max(1, ceil(n / 10));
  end
  if ~(isscalar(opts.progress_every) && isnumeric(opts.progress_every) && ...
       isfinite(opts.progress_every) && opts.progress_every == floor(opts.progress_every) && ...
       opts.progress_every >= 1)
    error('assemble_from_apply:InvalidProgressEvery: opts.progress_every must be an integer >= 1.');
  end

  if ~isfield(opts, 'force_full') || isempty(opts.force_full)
    opts.force_full = true;
  end
  if ~islogical(opts.force_full) || ~isscalar(opts.force_full)
    error('assemble_from_apply:InvalidForceFull: opts.force_full must be a logical scalar.');
  end

  % ----------------------------
  % Preallocation
  % ----------------------------
  % We store a dense matrix because later spectral computations (chol/eig) are
  % performed on dense matrices for small n.
  Mat = zeros(n, n);

  % Reuse one basis vector to avoid repeated allocations.
  e = zeros(n, 1);

  % ----------------------------
  % Column-by-column assembly
  % ----------------------------
  for j = 1:n
    % Form the j-th unit vector e_j.
    % (We reuse 'e' to reduce allocations; Octave uses copy-on-write semantics.)
    e(j) = 1;

    % Apply the matrix-free operator.
    y = applyFun(e);

    % Restore e for the next iteration.
    e(j) = 0;

    % Enforce column-vector shape and correct length.
    if ~isvector(y)
      error('assemble_from_apply:InvalidApplyOutput: applyFun returned a non-vector output at j=%d.', j);
    end
    y = y(:);

    if numel(y) ~= n
      error('assemble_from_apply:InvalidApplyOutputSize: applyFun output length %d != n=%d at j=%d.', numel(y), n, j);
    end

    % Ensure numeric type; allow sparse outputs but store dense if requested.
    if ~isnumeric(y) || ~isreal(y)
      error('assemble_from_apply:InvalidApplyOutputType: applyFun must return a real numeric vector (j=%d).', j);
    end
    if opts.force_full
      y = full(y);
    end

    % Store as the j-th column of the explicit matrix.
    Mat(:, j) = y;

    % Optional progress printing (off by default).
    if opts.verbose
      if (j == 1) || (j == n) || (mod(j, opts.progress_every) == 0)
        fprintf('[assemble_from_apply] assembled column %d / %d\n', j, n);
      end
    end
  end
end