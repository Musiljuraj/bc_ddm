% ============================================================
% File: main/ch6_define_cases.m
% ============================================================
function cases = ch6_define_cases()
%CH6_DEFINE_CASES  Define the Chapter 6 (sequential) experiment set (Set A + Set B).
%
% Returns:
%   cases : struct array, each entry has fields
%           - n      : mesh size parameter (NxN in your current setup)
%           - nSubX  : number of subdomains in x-direction
%           - nSubY  : number of subdomains in y-direction
%           - seed   : RNG seed (for reproducibility)
%
% Structure of the case set:
%   Set A (full spectra intended; moderate sizes):
%     - refinement at fixed decomposition:
%         (16,2x2) -> (24,2x2) -> (32,2x2) -> (48,2x2)
%     - decomposition effect at fixed mesh sizes:
%         (32,2x2) -> (32,4x4)
%         (48,2x2) -> (48,4x4)
%
%   Set B (weak-scaling-like hint; spectra may be skipped if too large):
%     - keep local mesh per subdomain roughly constant:
%         (32,4x4): 32/4 = 8 elements per subdomain side
%         (64,8x8): 64/8 = 8 elements per subdomain side
%
% Notes:
% - All n values are divisible by nSubX and nSubY (required by build_problem_data).
% - Keep seed fixed unless you explicitly want RHS variability (not needed here).

  seed = 1;

  cases = struct([]);

  k = 0;

  % ============================================================
  % Set A: full spectra intended
  % ============================================================

  % A1: baseline
  k = k + 1;
  cases(k).n     = 16;
  cases(k).nSubX = 2;
  cases(k).nSubY = 2;
  cases(k).seed  = seed;

  % A2: refinement (2x2)
  k = k + 1;
  cases(k).n     = 24;
  cases(k).nSubX = 2;
  cases(k).nSubY = 2;
  cases(k).seed  = seed;

  % A3: refinement (2x2)
  k = k + 1;
  cases(k).n     = 32;
  cases(k).nSubX = 2;
  cases(k).nSubY = 2;
  cases(k).seed  = seed;

  % A4: more subdomains at fixed mesh (32)
  k = k + 1;
  cases(k).n     = 32;
  cases(k).nSubX = 4;
  cases(k).nSubY = 4;
  cases(k).seed  = seed;

  % A5: refinement (2x2)
  k = k + 1;
  cases(k).n     = 48;
  cases(k).nSubX = 2;
  cases(k).nSubY = 2;
  cases(k).seed  = seed;

  % A6: more subdomains at fixed mesh (48)
  k = k + 1;
  cases(k).n     = 48;
  cases(k).nSubX = 4;
  cases(k).nSubY = 4;
  cases(k).seed  = seed;

  % ============================================================
  % Set B: weak-scaling-like hint (spectra may be skipped)
  % ============================================================

  % B1: weak-scaling-like larger case
  k = k + 1;
  cases(k).n     = 64;
  cases(k).nSubX = 8;
  cases(k).nSubY = 8;
  cases(k).seed  = seed;
end