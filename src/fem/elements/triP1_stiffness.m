function [Ke, area, gradphi] = triP1_stiffness(xy)
%TRIP1_STIFFNESS  Local stiffness matrix for Poisson problem on one triangular P1 element.
%
% Link to thesis:
%   Chapter 2.3.2 (Element stiffness matrix) and the element-wise decomposition
%   of K_ij = ∫_Ω ∇φ_i·∇φ_j dx.
%
% Theoretical role:
%   For each triangle T_e with local P1 basis functions {ϕ_1,ϕ_2,ϕ_3},
%     K^(e)_{ij} = ∫_{T_e} ∇ϕ_i · ∇ϕ_j dx.
%   This is the standard local (3x3) stiffness matrix for the bilinear form
%       a(u,v) = ∫_Ω  ∇(u)∇(v) dx
%   corresponding to the Poisson operator -Delta.
%   Since ϕ_i is affine on T_e, ∇ϕ_i is constant on T_e, so the integral reduces
%   to area(T_e) times the dot product of constant gradients.
%
%       Ke(i,j) = ( ∇φ_i·∇φ_j ) * |T_e|
%
%  Inputs:
%   xy      : 3×2 matrix -  coordinates of the triangle vertices.
%  Outputs:
%   Ke      : 3×3 local stiffness matrix for T_e.
%   area    : |T_e|, triangle area (positive)
%   gradphi : gradients of local basis functions (constant on the element).
%
% Implementation outline:
%   1) Compute the triangle area from the vertex geometry.
%   2) Construct the constant gradients of the local barycentric P1 basis.
%   3) Form Ke by area-scaled inner products of these gradients.
%   This routine is purely local and independent of global numbering.

  % -------------------------------------------------------------------------
  % 1) INPUT CHECKS
  % -------------------------------------------------------------------------
  if nargin ~= 1
    error('triP1_stiffness: expected exactly one input argument xy (3x2).');
  end

  % % CHANGED: stronger type/shape validation (numeric, real, finite, non-logical, 3x2)
  % Rationale:
  %   - Reject char arrays that look like matrices (ASCII codes).
  %   - Reject logical arrays (true/false treated as 1/0).
  %   - Reject NaN/Inf (otherwise degeneracy check can be bypassed).
  %   - Reject complex coordinates (not meaningful for geometry here).
  if ~ismatrix(xy) || any(size(xy) ~= [3, 2]) || ...
     ~isnumeric(xy) || islogical(xy) || ~isreal(xy) || any(~isfinite(xy(:)))
    error('triP1_stiffness:invalidInput', ...
          'xy must be a 3x2 real, finite, non-logical numeric matrix.');
  end

  % % ADDED: normalize to double for consistent downstream arithmetic
  xy = double(xy);

  % -------------------------------------------------------------------------
  % 2) READ VERTEX COORDINATES
  % -------------------------------------------------------------------------
  x1 = xy(1,1); y1 = xy(1,2);
  x2 = xy(2,1); y2 = xy(2,2);
  x3 = xy(3,1); y3 = xy(3,2);

  % -------------------------------------------------------------------------
  % 3) COMPUTE area of a parallelogram formed by two sides of an element (triangle)
  % -------------------------------------------------------------------------
  % Geometric meaning:
  %   Let v1 = [x2-x1; y2-y1], v2 = [x3-x1; y3-y1].
  %   det([v1 v2]) is the signed area of the parallelogram spanned by v1,v2.
  %   The element (triangle) area is half of that area.
  paralAreaSigned = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);

  % Degenerate triangle check (area ~ 0 => element is invalid for FEM).
  % NOTE: threshold remains absolute to preserve existing behavior in your pipeline/tests.
  if abs(paralAreaSigned) < 1e-14
    error('triP1_stiffness:degenerate', 'triP1_stiffness: degenerate triangle (area ~ 0).');
  end

  % Positive element area:
  absParalArea = abs(paralAreaSigned);     % % ADDED: reuse abs value
  area = 0.5 * absParalArea;

  % -------------------------------------------------------------------------
  % 4) COMPUTE b_i, c_i COEFFICIENTS (from nodal interpolation geometry)
  % -------------------------------------------------------------------------
  % For P1 basis functions on a triangle:
  %   phi_i(x,y) = a_i + b_i*x + c_i*y,
  % and the gradients are:
  %   grad(phi_i) = [b_i; c_i] / paralAreaSigned
  %
  % The b_i, c_i below are purely coordinate differences between the vertices.
  % They come from solving the interpolation conditions phi_i(x_j,y_j)=delta_ij for i = 1,2,3.
  b = [ y2 - y3;
        y3 - y1;
        y1 - y2 ];
  c = [ x3 - x2;
        x1 - x3;
        x2 - x1 ];

  % -------------------------------------------------------------------------
  % 5) GRADIENTS OF THE LOCAL BASIS FUNCTIONS
  % -------------------------------------------------------------------------
  % gradphi is 2x3:
  %   gradphi(:,i) is the constant gradient vector of basis function phi_i on K.
  %
  % Using signed paralAreaSigned keeps orientation consistent with vertex ordering.
  gradphi = [b.'; c.'] / paralAreaSigned;

  % -------------------------------------------------------------------------
  % 6) LOCAL STIFFNESS MATRIX
  % -------------------------------------------------------------------------
  % % CHANGED: compute Ke directly using an algebraically equivalent but slightly
  % more robust expression (avoids squaring paralAreaSigned).
  %
  % Starting from:
  %   Ke = area * (b*b' + c*c') / (paralAreaSigned^2)
  % and area = 0.5*abs(paralAreaSigned), we get:
  %   Ke = (b*b' + c*c') / (2*abs(paralAreaSigned)).
  Ke = (b*b.' + c*c.') / (2 * absParalArea);

end