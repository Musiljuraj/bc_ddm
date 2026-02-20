function sub = extract_subdomain_blocks(sub)
%EXTRACT_SUBDOMAIN_BLOCKS  Extract interior/interface blocks of local subdomain systems.
%
% Link to thesis:
%   Chapter 3.2.1 (block decomposition), equations (3.14)–(3.16).
%
% This routine assumes that each subdomain struct contains:
%   - K : local stiffness matrix on all local free DOFs
%   - f : local load vector
%   - dofs_I, dofs_G : local index sets (interior / interface)
%
% Output:
%   sub(i) is enriched with:
%     K_II, K_Ig, K_gI, K_gg
%     f_I,  f_g
%
% Notation:
%   "g" corresponds to Γ (interface) in the thesis (using g to avoid Unicode).

  if nargin ~= 1
    error('extract_subdomain_blocks: expected input (sub).');
  end

  nSub = numel(sub); 

  for i = 1:nSub
    if ~isfield(sub(i),'K') || ~isfield(sub(i),'f')
      error('extract_subdomain_blocks: sub(%d) must contain fields K and f.', i);
    end
    if ~isfield(sub(i),'dofs_I') || ~isfield(sub(i),'dofs_G')
      error('extract_subdomain_blocks: sub(%d) must contain fields dofs_I and dofs_G (run identify_interface_dofs first).', i);
    end

    I = sub(i).dofs_I(:); %interior nodes (in local, subdomain indexing, left-to-right, bottom-to-top)
    G = sub(i).dofs_G(:); %interface nodes (in local, subdomain indexing)

    K = sub(i).K; 
    f = sub(i).f;

    sub(i).K_II = K(I, I); %extracting block matrices and block RHS
    sub(i).K_Ig = K(I, G);
    sub(i).K_gI = K(G, I);
    sub(i).K_gg = K(G, G);

    sub(i).f_I  = f(I);
    sub(i).f_g  = f(G);

    sub(i).nI = numel(I);
    sub(i).nG = numel(G);
  end
end
