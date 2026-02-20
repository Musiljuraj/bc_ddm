function zC = apply_coarse_hook(data, r)
%APPLY_COARSE_HOOK  Placeholder for coarse (primal) correction application.
%
%   zC = apply_coarse_hook(data, r)
%
% Chapter 4.2 common interface ONLY:
%   The actual coarse correction is implemented in Chapters 4.3/4.4.
%   This placeholder fixes the signature early and is a safe no-op.
%
% Behavior:
%   returns zeros(size(r))

  %#ok<*INUSD>
  if nargin ~= 2
    error('apply_coarse_hook: expected inputs (data, r).');
  end
  zC = zeros(size(r));
end