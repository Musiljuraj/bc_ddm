% ============================================================
% File: src/feti_dp/solve/reconstruct_fetidp_solution.m
% ============================================================
function [w_c, w_d, u_free, diag] = reconstruct_fetidp_solution(lambda, data)
%RECONSTRUCT_FETIDP_SOLUTION  Recover interface traces and full free-DOF solution from lambda.
%
% Outputs:
%   w_c   : coarse vector u_c (assembled primal dofs, size nC)
%   w_d   : packed Delta vector (size nDeltaProd)
%   u_free: global solution in FREE-DOF numbering (size nFree x 1)
%   diag  : diagnostics (constraint residual norm, etc.)

  lambda = lambda(:);
  sub = data.sub;
  nSub = numel(sub);

  % Build gD packed and gC (data.hC) assembled into coarse ids using hc already computed.
  gD = zeros(data.nDeltaProd,1);
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    gi = sub(i).g(:);
    gD(rng) = gi(sub(i).idx_d(:));
  end

  % C*lambda contribution to coarse rhs:
  % Compute y = Sdd^{-1} (Bd^T lambda) and assemble -Scd*y into coarse rhs
  t = data.BdT * lambda;
  y = zeros(size(t));
  y_local = cell(nSub,1);
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      y_local{i} = zeros(0,1);
      continue;
    end
    ti = t(rng);
    R = data.Sdd_R{i};
    yi = R \ (R' \ ti);
    y(rng) = yi;
    y_local{i} = yi;
  end

  rhs_c = data.hc;
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      continue;
    end
    yi = y_local{i};
    if isempty(yi)
      continue;
    end
    rhs_c(c_ids) = rhs_c(c_ids) + data.Scd{i} * yi; % sign consistent with -Bd^T lambda inside delta equation
  end

  % Solve u_c
  if data.primal.nC > 0
    Rcc = data.Kcc_R;
    u_c = Rcc \ (Rcc' \ rhs_c);
  else
    u_c = zeros(0,1);
  end

  % Local primal traces wci
  wc_local = cell(nSub,1);
  for i = 1:nSub
    c_ids = sub(i).c_ids(:);
    if isempty(c_ids)
      wc_local{i} = zeros(0,1);
    else
      wc_local{i} = u_c(c_ids);
    end
  end

  % Compute delta trace:
  % wD = Sdd^{-1}( gD - Sdc*w_c - Bd^T lambda )
  rhsD = gD - t;
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    wci = wc_local{i};
    if ~isempty(wci)
      rhsD(rng) = rhsD(rng) - data.Sdc{i} * wci;
    end
  end

  wD = zeros(size(rhsD));
  for i = 1:nSub
    rng = data.delta_range{i};
    if isempty(rng)
      continue;
    end
    R = data.Sdd_R{i};
    wD(rng) = R \ (R' \ rhsD(rng));
  end

  % Reconstruct local gamma vectors and interior u_I via back-substitution
  nFree = numel(data.free);
  u_free_accum = zeros(nFree,1);
  u_free_count = zeros(nFree,1);

  for i = 1:nSub
    % local gamma in original sub(i).gamma ordering
    nG = numel(sub(i).idx_c) + numel(sub(i).idx_d);
    w_gamma = zeros(nG,1);
    % fill c and d parts
    w_gamma(sub(i).idx_c(:)) = wc_local{i};
    if ~isempty(sub(i).idx_d)
      w_gamma(sub(i).idx_d(:)) = wD(data.delta_range{i});
    end

    % Recover interior for this subdomain (if any)
    if isfield(sub(i), 'R_II') && ~isempty(sub(i).R_II)
      tmp = sub(i).f_I - sub(i).K_Ig * w_gamma;
      uI = sub(i).R_II \ (sub(i).R_II' \ tmp);
    else
      uI = zeros(0,1);
    end

    % Scatter to global FREE-DOF vector by averaging duplicates
    % Interface free dofs for this subdomain:
    globG = sub(i).glob_G(:);
    u_free_accum(globG) = u_free_accum(globG) + w_gamma;
    u_free_count(globG) = u_free_count(globG) + 1;

    % Interior free dofs:
    if isfield(sub(i), 'glob_I') && ~isempty(sub(i).glob_I)
      globI = sub(i).glob_I(:);
      u_free_accum(globI) = u_free_accum(globI) + uI;
      u_free_count(globI) = u_free_count(globI) + 1;
    end
  end

  % Final averaging for duplicates
  mask = (u_free_count > 0);
  u_free = zeros(nFree,1);
  u_free(mask) = u_free_accum(mask) ./ u_free_count(mask);

  % Outputs
  w_c = u_c;
  w_d = wD;

  diag = struct();
  diag.constraint_norm = norm(full(data.Bd * wD), 2);
end