function y = apply_blockdiag_S(sub, w)
%APPLY_BLOCKDIAG_S  Apply the block-diagonal interface operator S = diag(S_i) in the product space.
%
% Link to thesis:
%   Chapter 3.4.1, equation (3.54).
%
% After interior elimination, each subdomain contributes a local Schur complement S^(i)
% acting on its local interface DOFs. In the product (duplicated) space, the global
% operator is block-diagonal:
%     S = diag(S^(1), ..., S^(N)).
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

  if ~isstruct(sub)                                                % % ADDED
    error('apply_blockdiag_S: sub must be a struct array.');        % % ADDED
  end                                                              % % ADDED

  if ~isnumeric(w)                                                 % % ADDED
    error('apply_blockdiag_S: w must be numeric.');                 % % ADDED
  end                                                              % % ADDED
  if ~isreal(w)                                                    % % ADDED
    error('apply_blockdiag_S: w must be real-valued.');             % % ADDED
  end                                                              % % ADDED
  if any(~isfinite(w(:)))                                          % % ADDED
    error('apply_blockdiag_S: w must not contain NaN/Inf.');        % % ADDED
  end                                                              % % ADDED

  if ~isvector(w) || size(w,2) ~= 1
    error('apply_blockdiag_S: w must be a column vector.');
  end

  nProd = numel(w);
  y = zeros(nProd, 1);

  used = false(nProd, 1);                                          % % ADDED: detect overlaps in prod_idx

  nSub = numel(sub);
  for i = 1:nSub
    if ~isfield(sub(i),'prod_idx')
      error('apply_blockdiag_S: sub(%d) missing prod_idx (run build_product_interface).', i);
    end
    idx = sub(i).prod_idx(:);

    if isempty(idx)
      continue;
    end

    if any(idx ~= fix(idx))                                        % % ADDED
      error('apply_blockdiag_S: sub(%d).prod_idx must be integer-valued.', i); % % ADDED
    end                                                            % % ADDED

    if numel(unique(idx)) ~= numel(idx)                            % % ADDED
      error('apply_blockdiag_S: sub(%d).prod_idx contains duplicates.', i);   % % ADDED
    end                                                            % % ADDED

    if any(idx < 1) || any(idx > nProd)
      error('apply_blockdiag_S: sub(%d).prod_idx contains out-of-range indices.', i);
    end

    if any(used(idx))                                              % % ADDED
      error('apply_blockdiag_S: overlapping prod_idx across subdomains detected at sub(%d).', i); % % ADDED
    end                                                            % % ADDED
    used(idx) = true;                                              % % ADDED

    % Optional contract check: local interface size should match local K_gg dimension. % % ADDED
    if isfield(sub(i),'K_gg')                                      % % ADDED
      nG = size(sub(i).K_gg, 1);                                   % % ADDED
      if size(sub(i).K_gg,2) ~= nG                                 % % ADDED
        error('apply_blockdiag_S: sub(%d).K_gg must be square.', i);% % ADDED
      end                                                          % % ADDED
      if numel(idx) ~= nG                                          % % ADDED
        error('apply_blockdiag_S: sub(%d) size mismatch: |prod_idx|=%d but size(K_gg,1)=%d.', i, numel(idx), nG); % % ADDED
      end                                                          % % ADDED
    end                                                            % % ADDED

    wi = w(idx);
    yi = apply_local_schur(sub(i), wi);

    if ~isvector(yi) || size(yi,2) ~= 1 || numel(yi) ~= numel(idx)  % % ADDED
      error('apply_blockdiag_S: sub(%d) returned yi with unexpected shape/length.', i); % % ADDED
    end                                                            % % ADDED

    y(idx) = yi;
  end
end