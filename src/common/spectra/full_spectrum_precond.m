% ============================================================
% File: src/common/spectra/full_spectrum_precond.m
% ============================================================
function spec = full_spectrum_precond(applyA, applyMinv, n, opts)
%FULL_SPECTRUM_PRECOND Compute the full spectrum of a preconditioned SPD operator.
% Thesis link: Chapter 6.3 (explicit spectral analysis from operator actions).
% The routine works with matrices assembled from callbacks and returns the
% eigenvalues and derived spectral indicators used in Chapter 6.
%
% Inputs:
%   applyA     : function handle, y = applyA(x), x,y length n.
%   applyMinv  : function handle, z = applyMinv(r), r,z length n.
%   n          : positive integer, dimension of the iteration space.
%   opts       : optional struct with fields:
%                - symmetrize      (default true)  : replace A and Minv by their
%                                                   symmetric parts to suppress
%                                                   roundoff asymmetry.
%                - store_matrices  (default false) : store A/Minv/K in output.
%                - verbose         (default false) : print progress messages.
%                - assemble_opts   (default struct): options passed to
%                                                   assemble_from_apply.
%
% Output (struct spec):
%   spec.eigvals   : sorted eigenvalues of the preconditioned operator.
%   spec.lmin      : min eigenvalue
%   spec.lmax      : max eigenvalue
%   spec.kappa     : lmax / lmin (Inf if lmin <= 0)
%   spec.symmA     : relative symmetry defect of assembled A
%   spec.symmMinv  : relative symmetry defect of assembled Minv
%   spec.symmK     : relative symmetry defect of K
%   spec.chol_ok   : logical, true if Cholesky succeeded
%   spec.chol_msg  : string, diagnostic if Cholesky failed
%   (optional) spec.A, spec.Minv, spec.K if opts.store_matrices = true
%
% Notes:
% - This routine is intended for small n only (explicit assembly + full eig).
% - If chol(Minv) fails, eigenvalues are not computed and eigvals is empty.

  % ----------------------------
  % Input validation
  % ----------------------------
  if nargin < 3 || nargin > 4
    error('full_spectrum_precond:InvalidNargin: expected inputs (applyA, applyMinv, n [, opts]).');
  end

  if ~isa(applyA, 'function_handle')
    error('full_spectrum_precond:InvalidApplyA: applyA must be a function handle.');
  end

  if ~isa(applyMinv, 'function_handle')
    error('full_spectrum_precond:InvalidApplyMinv: applyMinv must be a function handle.');
  end

  if ~(isscalar(n) && isnumeric(n) && isreal(n) && isfinite(n) && n == floor(n) && n >= 1)
    error('full_spectrum_precond:InvalidN: n must be a positive integer scalar.');
  end

  if nargin < 4 || isempty(opts)
    opts = struct();
  end

  if ~isstruct(opts)
    error('full_spectrum_precond:InvalidOpts: opts must be a struct (or []).');
  end

  % ----------------------------
  % Options (defaults)
  % ----------------------------
  if ~isfield(opts, 'symmetrize') || isempty(opts.symmetrize)
    opts.symmetrize = true;
  end
  if ~islogical(opts.symmetrize) || ~isscalar(opts.symmetrize)
    error('full_spectrum_precond:InvalidSymmetrize: opts.symmetrize must be a logical scalar.');
  end

  if ~isfield(opts, 'store_matrices') || isempty(opts.store_matrices)
    opts.store_matrices = false;
  end
  if ~islogical(opts.store_matrices) || ~isscalar(opts.store_matrices)
    error('full_spectrum_precond:InvalidStoreMatrices: opts.store_matrices must be a logical scalar.');
  end

  if ~isfield(opts, 'verbose') || isempty(opts.verbose)
    opts.verbose = false;
  end
  if ~islogical(opts.verbose) || ~isscalar(opts.verbose)
    error('full_spectrum_precond:InvalidVerbose: opts.verbose must be a logical scalar.');
  end

  if ~isfield(opts, 'assemble_opts') || isempty(opts.assemble_opts)
    opts.assemble_opts = struct();
  end
  if ~isstruct(opts.assemble_opts)
    error('full_spectrum_precond:InvalidAssembleOpts: opts.assemble_opts must be a struct (or []).');
  end

  % If verbose is enabled here, enable coarse progress in assembly as well.
  if opts.verbose && ~isfield(opts.assemble_opts, 'verbose')
    opts.assemble_opts.verbose = true;
  end

  % ----------------------------
  % Initialize output structure
  % ----------------------------
  spec = struct();
  spec.eigvals   = [];
  spec.lmin      = NaN;
  spec.lmax      = NaN;
  spec.kappa     = NaN;

  spec.symmA     = NaN;
  spec.symmMinv  = NaN;
  spec.symmK     = NaN;

  spec.chol_ok   = false;
  spec.chol_msg  = '';

  % ----------------------------
  % Assemble explicit matrices from matrix-free actions
  % ----------------------------
  if opts.verbose
    fprintf('[full_spectrum_precond] assembling A (n=%d)\n', n);
  end
  A = assemble_from_apply(applyA, n, opts.assemble_opts);

  if opts.verbose
    fprintf('[full_spectrum_precond] assembling Minv (n=%d)\n', n);
  end
  Minv = assemble_from_apply(applyMinv, n, opts.assemble_opts);

  % ----------------------------
  % Symmetry diagnostics (before optional symmetrization)
  % ----------------------------
  spec.symmA    = rel_symm_defect_(A);
  spec.symmMinv = rel_symm_defect_(Minv);

  % ----------------------------
  % Optional symmetrization (recommended)
  % ----------------------------
  % In exact arithmetic, A and Minv should be symmetric in our setting, but
  % numerical roundoff may introduce small asymmetry. For stable eigenvalue
  % computations we remove this artifact by taking the symmetric parts.
  if opts.symmetrize
    A    = 0.5 * (A + A');
    Minv = 0.5 * (Minv + Minv');
  end

  % ----------------------------
  % Cholesky factorization of Minv
  % ----------------------------
  % Octave chol(Minv) returns upper-triangular R such that R'*R = Minv.
  if opts.verbose
    fprintf('[full_spectrum_precond] chol(Minv)\n');
  end

  try
    R = chol(Minv);
    spec.chol_ok  = true;
    spec.chol_msg = '';
  catch err
    spec.chol_ok  = false;
    spec.chol_msg = err.message;

    if opts.verbose
      fprintf('[full_spectrum_precond] chol failed: %s\n', spec.chol_msg);
    end

    % Store matrices for debugging if requested.
    if opts.store_matrices
      spec.A    = A;
      spec.Minv = Minv;
    end
    return;
  end

  % ----------------------------
  % Build symmetric matrix K and compute eigenvalues
  % ----------------------------
  if opts.verbose
    fprintf('[full_spectrum_precond] forming K = R*A*R''\n');
  end

  K = R * A * R';
  spec.symmK = rel_symm_defect_(K);

  if opts.verbose
    fprintf('[full_spectrum_precond] eig(K)\n');
  end

  eigvals = eig(K);

  % Eigenvalues should be real; protect against tiny imaginary parts.
  eigvals = real(eigvals(:));

  % Sort ascending for consistent plotting/tables.
  eigvals = sort(eigvals, 'ascend');

  spec.eigvals = eigvals;

  % ----------------------------
  % Derived scalars: lmin, lmax, kappa
  % ----------------------------
  if ~isempty(eigvals)
    spec.lmin = eigvals(1);
    spec.lmax = eigvals(end);

    if spec.lmin > 0
      spec.kappa = spec.lmax / spec.lmin;
    else
      spec.kappa = Inf;
    end
  end

  % ----------------------------
  % Optional storage of matrices (debugging / reproducibility)
  % ----------------------------
  if opts.store_matrices
    spec.A    = A;
    spec.Minv = Minv;
    spec.K    = K;
  end
end


% ============================================================
% Local helper: relative symmetry defect
% ============================================================
function s = rel_symm_defect_(M)
%REL_SYMM_DEFECT_  Relative symmetry defect ||M-M^T||_F / ||M||_F.
  nrm = norm(M, 'fro');
  if nrm == 0
    s = 0;
  else
    s = norm(M - M', 'fro') / nrm;
  end
end