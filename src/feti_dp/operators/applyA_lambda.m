% ============================================================
% File: src/feti_dp/operators/applyA_lambda.m
% ============================================================
function y = applyA_lambda(x, data)
%APPLYA_LAMBDA Apply the FETI-DP multiplier-space operator in matrix-free form.
% Thesis link: Chapter 5.3 and Chapter 6.2 (FETI-DP operator seen by PCG).
% The routine evaluates `A_lambda x` using `B_d^T`, `\tilde S^{-1}`, and `B_d`
% without assembling the global multiplier-space matrix explicitly.
  % ----------------------------
  % Input validation
  % ----------------------------
  if nargin ~= 2                                          
    error('applyA_lambda:InvalidNumArgs', ...              
          'Expected exactly 2 input arguments.');          
  end                                                      

  % Validate x
  if ~isnumeric(x) || islogical(x)                         
    error('applyA_lambda:InvalidXType', ...                
          'x must be a non-logical numeric vector.');      
  end                                                      
  if ~isreal(x)                                            
    error('applyA_lambda:InvalidXReal', ...                
          'x must be real-valued.');                       
  end                                                      
  if ~isvector(x)                                          
    error('applyA_lambda:InvalidXShape', ...               
          'x must be a vector.');                          
  end                                                      

  x = x(:);                                               

  if any(~isfinite(x))                                     
    error('applyA_lambda:InvalidXFinite', ...              
          'x must contain only finite values.');           
  end                                                      

  % Validate data
  if ~isstruct(data)                                       
    error('applyA_lambda:InvalidDataType', ...             
          'data must be a struct.');                       
  end                                                      
  if ~isfield(data, 'BdT') || ~isfield(data, 'Bd')          
    error('applyA_lambda:MissingDataFields', ...           
          'data must contain fields BdT and Bd.');         
  end                                                      

  BdT = data.BdT;                                          
  Bd  = data.Bd;                                           

  if ~isnumeric(BdT) || islogical(BdT)                     
    error('applyA_lambda:InvalidBdTType', ...              
          'data.BdT must be a numeric matrix/operator.');  
  end                                                      
  if ~isnumeric(Bd) || islogical(Bd)                       
    error('applyA_lambda:InvalidBdType', ...               
          'data.Bd must be a numeric matrix/operator.');   
  end                                                      

  nLambda = numel(x);                                      
  if size(BdT, 2) ~= nLambda                               
    error('applyA_lambda:DimMismatchBdT', ...              
          'Size mismatch: BdT must have size(*,%d).', nLambda);  
  end                                                      
  nDelta = size(BdT, 1);                                   
  if size(Bd, 2) ~= nDelta                                 
    error('applyA_lambda:DimMismatchBd', ...               
          'Size mismatch: Bd must have size(*,%d).', nDelta);    
  end                                                      
  if size(Bd, 1) ~= nLambda                                
    error('applyA_lambda:DimMismatchOutput', ...           
          'Size mismatch: Bd must have %d rows to match x.', nLambda);  
  end                                                      

  if exist('solve_tildeS', 'file') ~= 2                    
    error('applyA_lambda:MissingDependency', ...          
          'Required function solve_tildeS is not on the path.'); 
  end                                                    

  % ----------------------------
  % Core computation
  % ----------------------------
  t = BdT * x;                                            
  out = solve_tildeS([], t, data);                        

  if ~isstruct(out) || ~isfield(out, 'wD')                
    error('applyA_lambda:InvalidSolveOutput', ...         
          'solve_tildeS must return a struct with field wD.'); 
  end                                                     
  wD = out.wD;                                         
  if ~isnumeric(wD) || islogical(wD)                      
    error('applyA_lambda:InvalidWDType', ...              
          'solve_tildeS output wD must be numeric.');     
  end                                                     
  if ~isvector(wD) || numel(wD) ~= nDelta                 
    error('applyA_lambda:InvalidWDShape', ...             
          'solve_tildeS output wD must be a vector of length %d.', nDelta);
  end                                                     
  wD = wD(:);                                             
  if ~isreal(wD) || any(~isfinite(wD))               
    error('applyA_lambda:InvalidWDValues', ...         
          'solve_tildeS output wD must be real and finite.'); 
  end                                                     %

  y = Bd * wD;                                           
end