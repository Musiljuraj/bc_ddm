function y = apply_blockdiag_S(sub, w)
%APPLY_BLOCKDIAG_S Apply the block-diagonal Schur operator in product space.
% Thesis link: Chapter 4.2.4 and 4.3 (global interface operator structure).
% The routine applies local Schur complements independently on all subdomains.
%
% Inputs:
%   sub : subdomain array with fields:
%         .prod_idx (product indices for this subdomain's interface DOFs)
%         and the local Schur data required by apply_local_schur (K blocks + R_II).
%   w   : product interface vector (length nProd).
%
% Output:
%   y   : product interface vector y = S*w.

  if nargin ~= 2
    error('apply_blockdiag_S: expected inputs (sub, w).');
  end

  if ~isstruct(sub)                                               
    error('apply_blockdiag_S: sub must be a struct array.');        
  end                                                              

  if ~isnumeric(w)                                                 
    error('apply_blockdiag_S: w must be numeric.');                 
  end                                                              
  if ~isreal(w)                                                  
    error('apply_blockdiag_S: w must be real-valued.');           
  end                                                             
  if any(~isfinite(w(:)))                                        
    error('apply_blockdiag_S: w must not contain NaN/Inf.');        
  end                                                             
  if ~isvector(w) || size(w,2) ~= 1
    error('apply_blockdiag_S: w must be a column vector.');
  end

  nProd = numel(w);
  y = zeros(nProd, 1);

  used = false(nProd, 1);                                          

  nSub = numel(sub);
  for i = 1:nSub
    if ~isfield(sub(i),'prod_idx')
      error('apply_blockdiag_S: sub(%d) missing prod_idx (run build_product_interface).', i);
    end
    idx = sub(i).prod_idx(:);

    if isempty(idx)
      continue;
    end

    if any(idx ~= fix(idx))                                        
      error('apply_blockdiag_S: sub(%d).prod_idx must be integer-valued.', i); 
    end                                                            

    if numel(unique(idx)) ~= numel(idx)                         
      error('apply_blockdiag_S: sub(%d).prod_idx contains duplicates.', i);   
    end                                                            

    if any(idx < 1) || any(idx > nProd)
      error('apply_blockdiag_S: sub(%d).prod_idx contains out-of-range indices.', i);
    end

    if any(used(idx))                                              
      error('apply_blockdiag_S: overlapping prod_idx across subdomains detected at sub(%d).', i); 
    end                                                            
    used(idx) = true;                                              

    % Optional contract check: local interface size should match local K_gg dimension. 
    if isfield(sub(i),'K_gg')                                      
      nG = size(sub(i).K_gg, 1);                                   
      if size(sub(i).K_gg,2) ~= nG                                 
        error('apply_blockdiag_S: sub(%d).K_gg must be square.', i);
      end                                                          
      if numel(idx) ~= nG                                          
        error('apply_blockdiag_S: sub(%d) size mismatch: |prod_idx|=%d but size(K_gg,1)=%d.', i, numel(idx), nG); 
      end                                                         
    end                                                           

    wi = w(idx);
    yi = apply_local_schur(sub(i), wi);

    if ~isvector(yi) || size(yi,2) ~= 1 || numel(yi) ~= numel(idx)  
      error('apply_blockdiag_S: sub(%d) returned yi with unexpected shape/length.', i);
    end                                                            

    y(idx) = yi;
  end
end