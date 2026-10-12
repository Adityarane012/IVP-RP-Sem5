%% Pooled real-data threshold calibration
% Run from the code/ folder in MATLAB, after extracting the clean* crops.
%
% Requires in the same folder:
%   pooled_calibrate.m, rom_map.m, os_cfar.m, ordstat_filter.m, pad_symmetric.m

clear; clc;
import rosed.*          % the shared library, see code\rosed.m

CROPDIR = '../data/crops';
PATTERN = '*clean*.tif';
SCALES  = [8 14 20];
P_FA    = 1e-3;

[taus_real, alphas_real] = pooled_calibrate(CROPDIR, PATTERN, SCALES, P_FA);

if ~exist('../results/tables', 'dir'), mkdir('../results/tables'); end
save('../results/tables/calibration_real_pooled.mat', ...
     'taus_real', 'alphas_real', 'SCALES', 'P_FA', 'CROPDIR', 'PATTERN');

fprintf('\nSaved to ../results/tables/calibration_real_pooled.mat\n');
