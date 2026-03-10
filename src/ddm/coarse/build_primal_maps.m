function [sub, primal] = build_primal_maps(sub, ddm, iface, prod, primal)
%BUILD_PRIMAL_MAPS Build maps related to primal DOFs and coarse indexing.
% Thesis link: Chapter 4.4.2–4.4.3 (coarse interface space and reduced constraints).
% The routine prepares the local/global bookkeeping needed by dual-primal solvers.
%
% Inputs:
%   sub   : subdomain array with fields:
%           .glob_G or .gamma_glob : global free DOFs on interface (local order)
%           .prod_idx             : product indices (from build_product_interface)
%   ddm   : global bookkeeping (needs .nFree)
%   iface : interface bookkeeping (needs .glob2hat and .nHat)
%   prod  : product bookkeeping (needs .nProd, .prod2hat)
%   primal: primal selection (from select_primal_dofs), needs .glob_c (and optionally .hat_c)
%
% Outputs:
%   sub : updated with per-subdomain sets:
%       .idx_c         : indices of primal DOFs in the local interface vector (1..nG)
%       .idx_d         : indices of delta DOFs in the local interface vector (1..nG)
%       .glob_c        : global free DOFs of primal interface nodes on this subdomain
%       .glob_d        : global free DOFs of delta interface nodes on this subdomain
%       .c_ids         : global primal ids (1..nC) corresponding to glob_c
%       .prod_idx_c    : product indices corresponding to primal interface DOFs
%       .prod_idx_d    : product indices corresponding to delta interface DOFs
%
%   primal : enriched with global maps:
%       .glob2c        : (nFree x 1) map global free DOF -> primal id (0 if not primal)
%       .hat2c         : (nHat x 1)  map assembled interface index -> primal id (0 if not primal)
%       .prod_is_c     : (nProd x 1) logical, true if product DOF corresponds to a primal physical DOF
%       .prod_is_d     : (nProd x 1) logical complement
%       .delta_idx     : list of product indices belonging to delta part
%       .prod2delta    : (nProd x 1) map product index -> delta index (0 if primal)
%       .nDeltaProd    : number of delta product DOFs
%       .nC            : number of primal DOFs
%       .hat_c         : (nC x 1) hat indices of primal DOFs (normalized/verified)

  % -----------------------------
  % Basic arity check (keep simple; Octave may throw before body on wrong arity)
  % -----------------------------
  if nargin ~= 5
    error('build_primal_maps: expected inputs (sub, ddm, iface, prod, primal).');
  end

  % -----------------------------
  % Validate required structs/fields
  % -----------------------------
  if ~isstruct(ddm) || ~isfield(ddm,'nFree')
    error('build_primal_maps: ddm must be a struct containing nFree.');
  end
  if ~isstruct(iface) || ~isfield(iface,'nHat') || ~isfield(iface,'glob2hat')
    error('build_primal_maps: iface must be a struct containing nHat and glob2hat.');
  end
  if ~isstruct(prod) || ~isfield(prod,'nProd') || ~isfield(prod,'prod2hat')
    error('build_primal_maps: prod must be a struct containing nProd and prod2hat.');
  end
  if ~isstruct(primal) || ~isfield(primal,'glob_c')
    error('build_primal_maps: primal must be a struct containing glob_c (run select_primal_dofs).');
  end

  % -----------------------------
  % Harden sizes and index maps (fail fast on NaN/Inf/complex/non-integer/out-of-range)
  % -----------------------------
  nFree = local_require_int_scalar_(ddm.nFree, 'ddm.nFree');           
  nHat  = local_require_int_scalar_(iface.nHat, 'iface.nHat');      
  nProd = local_require_int_scalar_(prod.nProd, 'prod.nProd');        

  % iface.glob2hat: should map global free dof -> hat index (0 if not on interface).
  glob2hat = iface.glob2hat(:);                                        
  local_require_int_vector_nonneg_(glob2hat, 'iface.glob2hat');       
  if numel(glob2hat) ~= nFree                                          
    error('build_primal_maps: iface.glob2hat must have length nFree.');
  end
  if any(glob2hat > nHat)                                               
    error('build_primal_maps: iface.glob2hat contains values > nHat.');
  end

  % prod.prod2hat: should map product dof -> hat dof (must be in 1..nHat for existing products).
  prod2hat = prod.prod2hat(:);                                         
  local_require_int_vector_pos_(prod2hat, 'prod.prod2hat');             
  if numel(prod2hat) ~= nProd                                           
    error('build_primal_maps: prod.prod2hat must have length nProd.');
  end
  if nProd > 0 && nHat == 0                                           
    error('build_primal_maps: nHat==0 but nProd>0 (inconsistent interface/product bookkeeping).');
  end
  if nProd > 0
    if any(prod2hat < 1) || any(prod2hat > nHat)                        
      error('build_primal_maps: prod.prod2hat must be in 1..nHat.');
    end
  end

  % -----------------------------
  % Validate and normalize primal selection
  % -----------------------------
  glob_c = primal.glob_c(:);                                           
  local_require_int_vector_pos_(glob_c, 'primal.glob_c');             
  if any(glob_c < 1) || any(glob_c > nFree)                           
    error('build_primal_maps: primal.glob_c must be within 1..nFree.');
  end
  if numel(unique(glob_c)) ~= numel(glob_c)                          
    error('build_primal_maps: primal.glob_c contains duplicates (invalid primal set).');
  end

  nC = numel(glob_c);
  primal.nC = nC;

  if nC > 0 && nHat == 0                                              
    error('build_primal_maps: non-empty primal.glob_c but iface.nHat==0.');
  end

  % Compute hat_c from glob_c using glob2hat; primal DOFs must live on the interface.
  hat_c_from_glob = glob2hat(glob_c);                                  
  if any(hat_c_from_glob == 0)                                      
    error('build_primal_maps: primal.glob_c contains DOFs not present on the interface (glob2hat==0).');
  end
  if any(hat_c_from_glob < 1) || any(hat_c_from_glob > nHat)           
    error('build_primal_maps: iface.glob2hat(primal.glob_c) out of 1..nHat.');
  end
  if numel(unique(hat_c_from_glob)) ~= nC                               
    error('build_primal_maps: iface.glob2hat(primal.glob_c) has duplicates (inconsistent interface mapping).');
  end

  % If primal.hat_c is provided, verify it matches the implied mapping.
  if isfield(primal,'hat_c') && ~isempty(primal.hat_c)
    hat_c = primal.hat_c(:);                                          
    local_require_int_vector_pos_(hat_c, 'primal.hat_c');               
    if numel(hat_c) ~= nC                                             
      error('build_primal_maps: primal.hat_c must have same length as primal.glob_c.');
    end
    if any(hat_c < 1) || any(hat_c > nHat)                              
      error('build_primal_maps: primal.hat_c must be within 1..nHat.');
    end
    if numel(unique(hat_c)) ~= nC                                      
      error('build_primal_maps: primal.hat_c contains duplicates.');
    end
    if any(hat_c ~= hat_c_from_glob)                                    
      error('build_primal_maps: primal.hat_c is inconsistent with iface.glob2hat(primal.glob_c).');
    end
  else
    hat_c = hat_c_from_glob;                                            
  end
  primal.hat_c = hat_c;                                             

  % -----------------------------
  % Build global maps glob2c and hat2c
  % -----------------------------
  glob2c = zeros(nFree, 1);
  if nC > 0
    glob2c(glob_c) = (1:nC).';
  end
  primal.glob2c = glob2c;

  hat2c = zeros(nHat, 1);
  if nC > 0
    hat2c(hat_c) = (1:nC).';
  end
  primal.hat2c = hat2c;

  % -----------------------------
  % Build product primal/delta masks globally from prod2hat (not from sub loops)
  % -----------------------------
  if nProd == 0
    prod_is_c = false(0,1);
  else
    prod_is_c = (hat2c(prod2hat) > 0);                                  
  end
  primal.prod_is_c = prod_is_c;
  primal.prod_is_d = ~prod_is_c;

  delta_idx = find(~prod_is_c);
  primal.delta_idx = delta_idx;
  primal.nDeltaProd = numel(delta_idx);

  prod2delta = zeros(nProd, 1);
  prod2delta(delta_idx) = (1:primal.nDeltaProd).';
  primal.prod2delta = prod2delta;

  % -----------------------------
  % Per-subdomain split + cross-consistency checks
  % -----------------------------
  if ~isempty(sub) && ~isstruct(sub)                                
    error('build_primal_maps: sub must be a struct array (or empty).');
  end

  nSub = numel(sub);
  for i = 1:nSub

    % Allow either gamma_glob or glob_G; standardize to gamma_glob
    if ~isfield(sub(i),'gamma_glob')                                    
      if isfield(sub(i),'glob_G')
        sub(i).gamma_glob = sub(i).glob_G(:);
      else
        error('build_primal_maps: sub(%d) must contain glob_G or gamma_glob.', i);
      end
    end
    if ~isfield(sub(i),'prod_idx')
      error('build_primal_maps: sub(%d) missing prod_idx (run build_product_interface).', i);
    end

    gG = sub(i).gamma_glob(:);
    pidx = sub(i).prod_idx(:);

    % Validate gamma_glob and prod_idx
    local_require_int_vector_pos_(gG, sprintf('sub(%d).gamma_glob', i));     
    if any(gG < 1) || any(gG > nFree)                                      
      error('build_primal_maps: sub(%d).gamma_glob must be within 1..nFree.', i);
    end

    if numel(pidx) ~= numel(gG)                                             
      error('build_primal_maps: sub(%d) prod_idx length mismatch with gamma_glob.', i);
    end
    local_require_int_vector_pos_(pidx, sprintf('sub(%d).prod_idx', i));     
    if any(pidx < 1) || any(pidx > nProd)                                   
      error('build_primal_maps: sub(%d).prod_idx must be within 1..nProd.', i);
    end

    nG = numel(gG);
    if nG == 0
      sub(i).idx_c = zeros(0,1);
      sub(i).idx_d = zeros(0,1);
      sub(i).glob_c = zeros(0,1);
      sub(i).glob_d = zeros(0,1);
      sub(i).c_ids  = zeros(0,1);
      sub(i).prod_idx_c = zeros(0,1);
      sub(i).prod_idx_d = zeros(0,1);
      continue;
    end

    % Cross-check that these local interface DOFs truly lie on the interface and match product->hat mapping.
    hat_from_glob = glob2hat(gG);                                        
    if any(hat_from_glob == 0)                                         
      error('build_primal_maps: sub(%d).gamma_glob contains non-interface DOFs (glob2hat==0).', i);
    end
    hat_from_prod = prod2hat(pidx);                                       
    if any(hat_from_prod ~= hat_from_glob)                               
      error('build_primal_maps: sub(%d) inconsistent mapping: prod2hat(prod_idx) ~= glob2hat(gamma_glob).', i);
    end

    % Identify primal vs delta in this local ordering
    c_ids_local = glob2c(gG);               % 0 if not primal
    is_c = (c_ids_local > 0);

    sub(i).idx_c = find(is_c);
    sub(i).idx_d = find(~is_c);

    sub(i).glob_c = gG(is_c);
    sub(i).glob_d = gG(~is_c);

    sub(i).c_ids  = c_ids_local(is_c);

    sub(i).prod_idx_c = pidx(is_c);
    sub(i).prod_idx_d = pidx(~is_c);

    % Optional consistency: local product indices should agree with global prod_is_c
    if any(prod_is_c(sub(i).prod_idx_c) ~= true) || any(prod_is_c(sub(i).prod_idx_d) ~= false)  
      error('build_primal_maps: sub(%d) local primal/delta split inconsistent with global prod_is_c.', i);
    end
  end
end

% ============================================================
% Local validation helpers (Octave-friendly, no message matching required in tests)
% ============================================================

function n = local_require_int_scalar_(x, name)
% % ADDED
  if ~(isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x))
    error('build_primal_maps: %s must be a real finite numeric scalar.', name);
  end
  if x ~= round(x)
    error('build_primal_maps: %s must be integer-valued.', name);
  end
  if x < 0
    error('build_primal_maps: %s must be nonnegative.', name);
  end
  n = x;
end

function local_require_int_vector_nonneg_(v, name)
% % ADDED
  if isempty(v)
    return;
  end
  if ~(isnumeric(v) && isreal(v))
    error('build_primal_maps: %s must be real numeric.', name);
  end
  if any(~isfinite(v(:)))
    error('build_primal_maps: %s must be finite.', name);
  end
  if any(v(:) ~= round(v(:)))
    error('build_primal_maps: %s must be integer-valued.', name);
  end
  if any(v(:) < 0)
    error('build_primal_maps: %s must be >= 0.', name);
  end
end

function local_require_int_vector_pos_(v, name)
% % ADDED
  if isempty(v)
    return;
  end
  if ~(isnumeric(v) && isreal(v))
    error('build_primal_maps: %s must be real numeric.', name);
  end
  if any(~isfinite(v(:)))
    error('build_primal_maps: %s must be finite.', name);
  end
  if any(v(:) ~= round(v(:)))
    error('build_primal_maps: %s must be integer-valued.', name);
  end
  if any(v(:) <= 0)
    error('build_primal_maps: %s must be positive (>=1).', name);
  end
end