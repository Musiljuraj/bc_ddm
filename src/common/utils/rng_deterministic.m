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
    error('rng_deterministic: expected exactly one input (seed).');
  end
  if ~(isscalar(seed) && isfinite(seed) && seed == floor(seed))
    error('rng_deterministic: seed must be a finite integer scalar.');
  end

  rand('seed', seed);
  randn('seed', seed);
end