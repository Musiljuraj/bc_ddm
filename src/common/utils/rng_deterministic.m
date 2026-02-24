function rng_deterministic(seed)
%RNG_DETERMINISTIC  Set deterministic RNG state for Octave rand and randn.
%
%   rng_deterministic(seed)
%
% Chapter 4.2 common tool: tests and experiments should be reproducible.
%
% Behavior:
%   rand('seed', seed);
%   randn('seed', seed);

  if nargin ~= 1
    error('rng_deterministic:InvalidNargin: expected exactly one input (seed).'); % CHANGED
  end
  if nargout ~= 0
    error('rng_deterministic:InvalidNargout: no outputs are returned.'); % ADDED
  end

  % Conservative hardening: reject non-numeric / logical / complex / non-finite. % ADDED
  if ~(isnumeric(seed) && isreal(seed) && isscalar(seed))
    error('rng_deterministic:InvalidSeed: seed must be a real numeric scalar.'); % CHANGED
  end
  if islogical(seed)
    error('rng_deterministic:InvalidSeed: logical seeds are not allowed; use an integer numeric type.'); % ADDED
  end

  s = double(seed); % ADDED (normalize for checks + seeding)

  if ~isfinite(s)
    error('rng_deterministic:InvalidSeed: seed must be finite.'); % CHANGED
  end
  if abs(s) > flintmax
    error('rng_deterministic:InvalidSeed: seed magnitude too large to be represented exactly as an integer.'); % ADDED
  end
  if s ~= fix(s)
    error('rng_deterministic:InvalidSeed: seed must be an integer-valued scalar.'); % CHANGED
  end

  rand('seed', s);   % CHANGED (seed using validated double)
  randn('seed', s);  % CHANGED
end