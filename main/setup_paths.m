% ============================================================
% File: main/setup_paths.m
% ============================================================
function setup_paths()
%SETUP_PATHS  Add project paths (src + tests) regardless of current working directory.
%
% This helper should be called at the top of every main/*.m entry point.

  here = fileparts(mfilename('fullpath'));     % .../main
  root = fileparts(here);                      % project root

  addpath(genpath(fullfile(root, 'src')));
  addpath(genpath(fullfile(root, 'tests')));
end