function primal = select_primal_dofs(p, ddm, iface)
%SELECT_PRIMAL_DOFS Select primal interface DOFs for the coarse space.
% Thesis link: Chapter 4.4.1–4.4.2 (primal variables and coarse continuity).
% In the structured setting used here, the routine selects the interface
% DOFs that remain globally identified.
%
% Inputs:
%   p    : Np x 2 node coordinates.                
%   ddm  : struct from build_subdomains_structured 
%          (needs nSubX,nSubY,dof2node).
%   iface: struct from identify_interface_dofs    
%          (needs glob and glob2hat).
%
% Output:
%   primal : struct with fields:
%       .glob_c     : global free DOFs selected as primal
%       .glob_delta : remaining interface global free DOFs
%       .hat_c      : assembled interface indices of primal DOFs
%       .hat_delta  : assembled interface indices of remaining interface DOFs
%       .nC         : number of primal DOFs
%       .tol        : tolerance used (for diagnostics)

  if nargin ~= 3
    error('select_primal_dofs: expected inputs (p, ddm, iface).');
  end

  if nargout > 1                                                  
    error('select_primal_dofs: too many output arguments.');       
  end

  % Fixed coordinate tolerance (no opts input).
  % CHANGED: scale-aware lower bound while preserving baseline 1e-12 on unit-square meshes.
  p_scale = max(1, max(abs(double(p(:)))));                       
  tol = max(1e-12, 50 * eps(p_scale));                             

  % CHANGED: strong input validation (fail fast vs silent ASCII/NaN propagation).
  validate_inputs_(p, ddm, iface, tol);                             

  nSubX = ddm.nSubX;
  nSubY = ddm.nSubY;

  % Lines that define subdomain boundaries (their cross-points are corners).
  xLines = (0:nSubX) / nSubX;
  yLines = (0:nSubY) / nSubY;

  iface_glob = iface.glob(:);
  nHat = numel(iface_glob);

  is_primal = false(nHat, 1);

  % For all dofs in w-hat find out whether they lay on a cross-point of subdomain boundary lines.
  for k = 1:nHat
    g = iface_glob(k);                 % global free DOF id
    node = ddm.dof2node(g);            % global node index (free dof -> node map)
    x = p(node,1);
    y = p(node,2);

    on_x = any(abs(x - xLines) < tol);
    on_y = any(abs(y - yLines) < tol);

    if on_x && on_y
      is_primal(k) = true;
    end
  end

  primal = struct();
  primal.glob_c = sort(iface_glob(is_primal));
  primal.glob_delta = sort(iface_glob(~is_primal));

  primal.hat_c = iface.glob2hat(primal.glob_c);
  primal.hat_delta = iface.glob2hat(primal.glob_delta);

  % validate glob2hat results are usable indices into iface.glob (catch zeros/out-of-range).
  if ~isempty(primal.hat_c)                                        
    if any(primal.hat_c < 1 | primal.hat_c > nHat)                 
      error('select_primal_dofs: iface.glob2hat returned out-of-range indices for glob_c.'); 
    end
  end
  if ~isempty(primal.hat_delta)                                   
    if any(primal.hat_delta < 1 | primal.hat_delta > nHat)         
      error('select_primal_dofs: iface.glob2hat returned out-of-range indices for glob_delta.');
    end
  end

  primal.nC = numel(primal.glob_c);

  % Store for convenience.
  primal.tol = tol;
end


% ======================= local helpers =======================

function validate_inputs_(p, ddm, iface, tol)
%VALIDATE_INPUTS_  Input contract checks for select_primal_dofs.


  % ---- p validation ----
  if ~isnumeric(p) || islogical(p)                             
    error('select_primal_dofs: p must be a numeric (non-logical) array.');
  end
  if ~isreal(p)                                                     
    error('select_primal_dofs: p must be real-valued.');            
  end
  if any(~isfinite(p(:)))                                           
    error('select_primal_dofs: p must be finite (no NaN/Inf).');    
  end
  if ndims(p) ~= 2 || size(p,2) ~= 2 || size(p,1) < 1               
    error('select_primal_dofs: p must be an Np x 2 coordinate array with Np>=1.'); 
  end
  if ~isscalar(tol) || ~isfinite(tol) || tol <= 0                  
    error('select_primal_dofs: internal error (tol must be positive finite scalar).'); 
  end

  % ---- ddm validation ----
  if ~isstruct(ddm)                                                 
    error('select_primal_dofs: ddm must be a struct.');             
  end
  if ~isfield(ddm,'nSubX') || ~isfield(ddm,'nSubY') || ~isfield(ddm,'dof2node')
    error('select_primal_dofs: ddm must contain nSubX, nSubY, dof2node.');
  end

  nSubX = ddm.nSubX;
  nSubY = ddm.nSubY;

  validate_pos_int_scalar_(nSubX, 'ddm.nSubX');                    
  validate_pos_int_scalar_(nSubY, 'ddm.nSubY');                     

  dof2node = ddm.dof2node;
  if ~isnumeric(dof2node) || islogical(dof2node) || ~isreal(dof2node)
    error('select_primal_dofs: ddm.dof2node must be a real numeric vector.'); 
  end
  if isempty(dof2node) || ~(isvector(dof2node))
    error('select_primal_dofs: ddm.dof2node must be a nonempty vector.');    
  end
  if any(~isfinite(dof2node(:)))
    error('select_primal_dofs: ddm.dof2node must be finite (no NaN/Inf).');  
  end
  if any(dof2node(:) ~= round(dof2node(:)))
    error('select_primal_dofs: ddm.dof2node must be integer-valued indices.'); 
  end
  nNodes = size(p,1);
  if any(dof2node(:) < 1 | dof2node(:) > nNodes)
    error('select_primal_dofs: ddm.dof2node indices must be in 1..size(p,1).'); 
  end

  % ---- iface validation ----
  if ~isstruct(iface)                                             
    error('select_primal_dofs: iface must be a struct.');         
  end
  if ~isfield(iface,'glob') || ~isfield(iface,'glob2hat')
    error('select_primal_dofs: iface must contain glob and glob2hat.');
  end

  glob = iface.glob;
  if ~isnumeric(glob) || islogical(glob) || ~isreal(glob)
    error('select_primal_dofs: iface.glob must be a real numeric vector.');   
  end
  if isempty(glob) || ~isvector(glob)
    error('select_primal_dofs: iface.glob must be a nonempty vector.');       
  end
  if any(~isfinite(glob(:)))
    error('select_primal_dofs: iface.glob must be finite (no NaN/Inf).');     
  end
  if any(glob(:) ~= round(glob(:)))
    error('select_primal_dofs: iface.glob must be integer-valued indices.');  
  end

  glob = glob(:);  % ensure column
  if numel(unique(glob)) ~= numel(glob)
    error('select_primal_dofs: iface.glob must not contain duplicates.');     
  end

  nFree = numel(dof2node);
  if any(glob < 1 | glob > nFree)
    error('select_primal_dofs: iface.glob entries must be in 1..numel(ddm.dof2node).'); 
  end

  glob2hat = iface.glob2hat;
  if ~isnumeric(glob2hat) || islogical(glob2hat) || ~isreal(glob2hat)
    error('select_primal_dofs: iface.glob2hat must be a real numeric array.'); 
  end
  if any(~isfinite(glob2hat(:)))
    error('select_primal_dofs: iface.glob2hat must be finite (no NaN/Inf).');  
  end
  if any(glob2hat(:) ~= round(glob2hat(:)))
    error('select_primal_dofs: iface.glob2hat must be integer-valued.');       
  end
  if numel(glob2hat) < max(glob)
    error('select_primal_dofs: iface.glob2hat must be indexable at max(iface.glob).'); 
  end

  % Mapping consistency: for every g in iface.glob, glob2hat(g) must be a unique hat index in 1..nHat.
  nHat = numel(glob);
  hat_of_glob = glob2hat(glob);                                   
  if any(hat_of_glob < 1 | hat_of_glob > nHat)
    error('select_primal_dofs: iface.glob2hat(iface.glob) must be in 1..numel(iface.glob).'); 
  end
  if numel(unique(hat_of_glob)) ~= nHat
    error('select_primal_dofs: iface.glob2hat mapping must be one-to-one on iface.glob.'); 
  end
  % Strongest consistency check with the current iface.glob ordering:
  if ~isequal(glob(hat_of_glob(:)), glob(:))
    error('select_primal_dofs: inconsistent iface.glob and iface.glob2hat (not inverse on iface.glob).'); 
  end

end

function validate_pos_int_scalar_(x, name)
%VALIDATE_POS_INT_SCALAR_  Require positive integer scalar.

  if ~isnumeric(x) || islogical(x) || ~isreal(x) || ~isscalar(x) || ~isfinite(x)
    error(sprintf('select_primal_dofs: %s must be a real finite numeric scalar.', name));
  end
  if x <= 0 || x ~= round(x)
    error(sprintf('select_primal_dofs: %s must be a positive integer scalar.', name));
  end
end