function fe = triP1_load(xy, f_handle)
%TRIP1_LOAD  Local load vector on one P1 triangular element (Poisson RHS).
%
% Link to thesis:
%   Chapter 2.3.1 (discrete RHS entries f_i = ∫_Ω φ_i f dx) and the element-wise
%   decomposition into triangle contributions.
%
% Theoretical role:
%   For one element T_e with local basis {ϕ_1,ϕ_2,ϕ_3}, we define the local element load vector
%     f^(e)_m = ∫_{T_e} ϕ_m(x) f(x) dx,  m=1,2,3,
%   and the global load vector is obtained by summing these local contributions
%   according to mesh connectivity (assembly).
%
% Quadrature model:
%   The element integral is approximated by a centroid one-point rule.
%   For P1 basis functions, ϕ_m(x_c)=1/3 at the centroid x_c, giving the compact
%   approximation f^(e) ≈ |T_e| f(x_c) (1/3)[1;1;1].
%
% Inputs:
%   xy       : 3×2 coordinates of triangle vertices. 
%   f_handle : function handle for the load (force) term f(x,y).
% Outputs:
%   fe       : 3×1 local load vector for element T_e.
%
% Implementation outline:
%   1) Compute area(T_e).
%   2) Evaluate f at the triangle centroid.
%   3) Form the centroid quadrature-based local vector fe.

  % -----------------------------
  % 1) Input validation
  % -----------------------------
  if nargin ~= 2
    error('triP1_load: expected inputs (xy, f_handle).');
  end
  if ~ismatrix(xy) || size(xy,1) ~= 3 || size(xy,2) ~= 2
    error('triP1_load: xy must be 3x2 (triangle vertex coordinates).');
  end
  if ~isa(f_handle, 'function_handle')
    error('triP1_load: f_handle must be a function handle, e.g. @(x,y) ...');
  end

  % -----------------------------
  % 2) Compute triangle area |T|
  % -----------------------------
  % The area of the triangle with vertices (x1,y1), (x2,y2), (x3,y3) is
  %   |T| = 0.5 * |det([x2-x1, x3-x1; y2-y1, y3-y1])|.
  %
  % This is the same geometric determinant logic you used in triP1_stiffness.
  x1 = xy(1,1); y1 = xy(1,2);
  x2 = xy(2,1); y2 = xy(2,2);
  x3 = xy(3,1); y3 = xy(3,2);
  
  % -------------------------------------------------------------------------
  % COMPUTE area of a parallelogram formed by two sides of an element (triangle)
  % -------------------------------------------------------------------------
  % Geometric meaning:
  %   Let v1 = [x2-x1; y2-y1], v2 = [x3-x1; y3-y1].
  %   det([v1 v2]) is the signed area of the parallelogram spanned by v1,v2.
  %   The element (triangle) area is half of that area.
  paralAreaSigned = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
  area = 0.5 * abs(paralAreaSigned);

  if area <= 0
    error('triP1_load: degenerate triangle (area <= 0).');
  end

  % -----------------------------
  % 3) Evaluate f at the centroid
  % -----------------------------
  % Triangle centroid:
  %   x_c = (x1+x2+x3)/3, y_c = (y1+y2+y3)/3
  xc = (x1 + x2 + x3) / 3;
  yc = (y1 + y2 + y3) / 3;

  % Evaluate the source term. This must return a scalar.
  fval = f_handle(xc, yc);
  if ~isscalar(fval) || ~isfinite(fval)
    error('triP1_load: f_handle must return a finite scalar.');
  end

  % -----------------------------
  % 4) 1-point quadrature formula
  % -----------------------------
  % Using φ_m(x_c) = 1/3 for P1 basis functions:
  %   fe(m) ≈ f(x_c) * area * (1/3)
  fe = (fval * area / 3) * ones(3,1);

end
