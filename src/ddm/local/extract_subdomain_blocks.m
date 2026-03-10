function sub = extract_subdomain_blocks(sub)
%EXTRACT_SUBDOMAIN_BLOCKS Split each local system into interior/interface blocks.
% Thesis link: Chapter 4.2.1 (block decomposition of local subdomain systems).
% The routine extracts `K_II`, `K_Ig`, `K_gI`, `K_gg` and the matching RHS parts.
%
% Output:
%   sub(i) is enriched with:
%     K_II, K_Ig, K_gI, K_gg
%     f_I,  f_g


  if nargin ~= 1
    error('extract_subdomain_blocks: expected input (sub).');
  end

  nSub = numel(sub);

  for i = 1:nSub
    % -----------------------------
    % % ADDED: required field checks (already present, kept)
    % -----------------------------
    if ~isfield(sub(i),'K') || ~isfield(sub(i),'f')
      error('extract_subdomain_blocks: sub(%d) must contain fields K and f.', i);
    end
    if ~isfield(sub(i),'dofs_I') || ~isfield(sub(i),'dofs_G')
      error('extract_subdomain_blocks: sub(%d) must contain fields dofs_I and dofs_G (run identify_interface_dofs first).', i);
    end

    % -----------------------------
    % % ADDED: validate K and f basic invariants
    % -----------------------------
    K = sub(i).K;
    f = sub(i).f;

    if ~isnumeric(K) || ~isreal(K)
      error('extract_subdomain_blocks: sub(%d).K must be a numeric real matrix.', i);
    end
    if ndims(K) ~= 2 || size(K,1) ~= size(K,2)
      error('extract_subdomain_blocks: sub(%d).K must be a square 2D matrix.', i);
    end

    n = size(K,1);

    % % ADDED: finite check (sparse-friendly)
    if issparse(K)
      nz = nonzeros(K);
      if ~all(isfinite(nz))
        error('extract_subdomain_blocks: sub(%d).K contains NaN/Inf entries.', i);
      end
    else
      if ~all(isfinite(K(:)))
        error('extract_subdomain_blocks: sub(%d).K contains NaN/Inf entries.', i);
      end
    end

    if ~isnumeric(f) || ~isreal(f) || ~isvector(f)
      error('extract_subdomain_blocks: sub(%d).f must be a numeric real vector.', i);
    end
    f = f(:);  % % ADDED: enforce column vector internally
    if numel(f) ~= n
      error('extract_subdomain_blocks: sub(%d).f must have length size(K,1).', i);
    end
    if ~all(isfinite(f))
      error('extract_subdomain_blocks: sub(%d).f contains NaN/Inf entries.', i);
    end

    % -----------------------------
    % % ADDED: validate dofs_I / dofs_G (reject char/logical; integer/range; duplicates; overlap)
    % -----------------------------
    Iraw = sub(i).dofs_I;
    Graw = sub(i).dofs_G;

    if ischar(Iraw) || islogical(Iraw)
      error('extract_subdomain_blocks: sub(%d).dofs_I must be numeric (reject char/logical).', i);
    end
    if ischar(Graw) || islogical(Graw)
      error('extract_subdomain_blocks: sub(%d).dofs_G must be numeric (reject char/logical).', i);
    end

    if ~isnumeric(Iraw) || ~isreal(Iraw) || ~(isvector(Iraw) || isempty(Iraw))
      error('extract_subdomain_blocks: sub(%d).dofs_I must be a numeric real vector (or empty).', i);
    end
    if ~isnumeric(Graw) || ~isreal(Graw) || ~(isvector(Graw) || isempty(Graw))
      error('extract_subdomain_blocks: sub(%d).dofs_G must be a numeric real vector (or empty).', i);
    end

    I = Iraw(:);  % % CHANGED: validate before flattening
    G = Graw(:);

    if ~all(isfinite(I)) || ~all(isfinite(G))
      error('extract_subdomain_blocks: sub(%d).dofs_I/dofs_G must be finite (no NaN/Inf).', i);
    end

    % % ADDED: integer-valued indices
    if ~isempty(I) && any(mod(I,1) ~= 0)
      error('extract_subdomain_blocks: sub(%d).dofs_I must be integer-valued.', i);
    end
    if ~isempty(G) && any(mod(G,1) ~= 0)
      error('extract_subdomain_blocks: sub(%d).dofs_G must be integer-valued.', i);
    end

    % % ADDED: index range
    if ~isempty(I) && (any(I < 1) || any(I > n))
      error('extract_subdomain_blocks: sub(%d).dofs_I indices must satisfy 1 <= idx <= size(K,1).', i);
    end
    if ~isempty(G) && (any(G < 1) || any(G > n))
      error('extract_subdomain_blocks: sub(%d).dofs_G indices must satisfy 1 <= idx <= size(K,1).', i);
    end

    % % ADDED: duplicates within each set (reject)
    if numel(unique(I)) ~= numel(I)
      error('extract_subdomain_blocks: sub(%d).dofs_I must not contain duplicates.', i);
    end
    if numel(unique(G)) ~= numel(G)
      error('extract_subdomain_blocks: sub(%d).dofs_G must not contain duplicates.', i);
    end

    % % ADDED: disjointness between interior and interface (reject overlap)
    if ~isempty(I) && ~isempty(G) && ~isempty(intersect(I, G))
      error('extract_subdomain_blocks: sub(%d).dofs_I and dofs_G must be disjoint.', i);
    end

    % -----------------------------
    % Block extraction (unchanged mathematically)
    % -----------------------------
    sub(i).K_II = K(I, I);
    sub(i).K_Ig = K(I, G);
    sub(i).K_gI = K(G, I);
    sub(i).K_gg = K(G, G);

    sub(i).f_I  = f(I);
    sub(i).f_g  = f(G);

    sub(i).nI = numel(I);
    sub(i).nG = numel(G);
  end
end