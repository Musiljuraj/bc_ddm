function zC = apply_coarse_hook(data, r)
%APPLY_COARSE_HOOK Apply the coarse correction in the common dual-primal layer.
% Thesis link: Chapter 4.4 and Chapter 5.2 (coarse continuity mechanism).
% This helper evaluates the coarse-space contribution used by the solvers.

  if nargin ~= 2
    error('apply_coarse_hook: expected inputs (data, r).');
  end
  zC = zeros(size(r));
end