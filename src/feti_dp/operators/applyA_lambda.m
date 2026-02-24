% ============================================================
% File: src/feti_dp/operators/applyA_lambda.m
% ============================================================
function y = applyA_lambda(x, data)
%APPLYA_LAMBDA  Matrix-free application of the reduced SPD operator A_lambda.
%
% y = A_lambda * x
%
% Operational form:
%   1) t = Bd^T x
%   2) wD = [ \tilde S^{-1} (0, t) ]_Delta
%   3) y = Bd * wD

  % ----------------------------
  % Input validation
  % ----------------------------
  if nargin ~= 2                                         % ADDED
    error('applyA_lambda:InvalidNumArgs', ...             % ADDED
          'Expected exactly 2 input arguments.');         % ADDED
  end                                                     % ADDED

  % Validate x
  if ~isnumeric(x) || islogical(x)                        % ADDED
    error('applyA_lambda:InvalidXType', ...               % ADDED
          'x must be a non-logical numeric vector.');     % ADDED
  end                                                     % ADDED
  if ~isreal(x)                                           % ADDED
    error('applyA_lambda:InvalidXReal', ...               % ADDED
          'x must be real-valued.');                      % ADDED
  end                                                     % ADDED
  if ~isvector(x)                                         % ADDED
    error('applyA_lambda:InvalidXShape', ...              % ADDED
          'x must be a vector.');                         % ADDED
  end                                                     % ADDED

  x = x(:);                                               % CHANGED (kept behavior; now after checks)

  if any(~isfinite(x))                                    % ADDED
    error('applyA_lambda:InvalidXFinite', ...             % ADDED
          'x must contain only finite values.');          % ADDED
  end                                                     % ADDED

  % Validate data
  if ~isstruct(data)                                      % ADDED
    error('applyA_lambda:InvalidDataType', ...            % ADDED
          'data must be a struct.');                      % ADDED
  end                                                     % ADDED
  if ~isfield(data, 'BdT') || ~isfield(data, 'Bd')         % ADDED
    error('applyA_lambda:MissingDataFields', ...          % ADDED
          'data must contain fields BdT and Bd.');        % ADDED
  end                                                     % ADDED

  BdT = data.BdT;                                         % ADDED
  Bd  = data.Bd;                                          % ADDED

  if ~isnumeric(BdT) || islogical(BdT)                    % ADDED
    error('applyA_lambda:InvalidBdTType', ...             % ADDED
          'data.BdT must be a numeric matrix/operator.'); % ADDED
  end                                                     % ADDED
  if ~isnumeric(Bd) || islogical(Bd)                      % ADDED
    error('applyA_lambda:InvalidBdType', ...              % ADDED
          'data.Bd must be a numeric matrix/operator.');  % ADDED
  end                                                     % ADDED

  nLambda = numel(x);                                     % ADDED
  if size(BdT, 2) ~= nLambda                              % ADDED
    error('applyA_lambda:DimMismatchBdT', ...             % ADDED
          'Size mismatch: BdT must have size(*,%d).', nLambda); % ADDED
  end                                                     % ADDED
  nDelta = size(BdT, 1);                                  % ADDED
  if size(Bd, 2) ~= nDelta                                % ADDED
    error('applyA_lambda:DimMismatchBd', ...              % ADDED
          'Size mismatch: Bd must have size(*,%d).', nDelta);   % ADDED
  end                                                     % ADDED
  if size(Bd, 1) ~= nLambda                               % ADDED
    error('applyA_lambda:DimMismatchOutput', ...          % ADDED
          'Size mismatch: Bd must have %d rows to match x.', nLambda); % ADDED
  end                                                     % ADDED

  if exist('solve_tildeS', 'file') ~= 2                   % ADDED
    error('applyA_lambda:MissingDependency', ...          % ADDED
          'Required function solve_tildeS is not on the path.'); % ADDED
  end                                                     % ADDED

  % ----------------------------
  % Core computation
  % ----------------------------
  t = BdT * x;                                            % CHANGED (use local BdT)
  out = solve_tildeS([], t, data);                        % unchanged call convention

  if ~isstruct(out) || ~isfield(out, 'wD')                % ADDED
    error('applyA_lambda:InvalidSolveOutput', ...         % ADDED
          'solve_tildeS must return a struct with field wD.'); % ADDED
  end                                                     % ADDED
  wD = out.wD;                                            % ADDED
  if ~isnumeric(wD) || islogical(wD)                      % ADDED
    error('applyA_lambda:InvalidWDType', ...              % ADDED
          'solve_tildeS output wD must be numeric.');     % ADDED
  end                                                     % ADDED
  if ~isvector(wD) || numel(wD) ~= nDelta                 % ADDED
    error('applyA_lambda:InvalidWDShape', ...             % ADDED
          'solve_tildeS output wD must be a vector of length %d.', nDelta); % ADDED
  end                                                     % ADDED
  wD = wD(:);                                             % ADDED
  if ~isreal(wD) || any(~isfinite(wD))                    % ADDED
    error('applyA_lambda:InvalidWDValues', ...            % ADDED
          'solve_tildeS output wD must be real and finite.'); % ADDED
  end                                                     % ADDED

  y = Bd * wD;                                            % CHANGED (use local Bd, validated wD)
end