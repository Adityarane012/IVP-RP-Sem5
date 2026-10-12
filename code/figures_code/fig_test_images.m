%% ========================================================================
%  fig_test_images.m  —  All five simulated test images in one figure
%
%  WHY THIS EXISTS
%   The paper used two figures for the test images, one for A, B and C and
%   one for D and E. They are the same kind of picture and were merged into
%   a single row to save a page.
%
%  WHAT IT USES (nothing is recomputed, no new random numbers)
%   ../data/simulated/fig6abcd.mat      written by rosed_pipeline.m
%   ../data/simulated/fig6e_images.mat  written by fig6e_extended.m
%
%  OUTPUT  ../results/figures/paper/fig02_test_images.png
%
%  Run from code/ with code/figures_code on the path.
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

SIM = '../data/simulated';
OUT = '../results/figures/paper';
if ~exist(OUT, 'dir'), mkdir(OUT); end

D = load(fullfile(SIM, 'fig6abcd.mat'));
E = load(fullfile(SIM, 'fig6e_images.mat'));
i10 = find(E.GEN_ORDER == 10, 1);
if isempty(i10), error('The 10-pixel separation is not in fig6e_images.mat'); end

imgs = {D.Pa, D.Pb, D.Pc, D.Pd, E.imgs{i10}};
ttl  = {'(a) A: clutter boundary', '(b) B: bright line', '(c) C: dark line', ...
        '(d) D: all three types', '(e) E: two boundaries, 10 px'};

f = paper_fig('Test images', [50 50 1750 420]);
for k = 1:numel(imgs)
    show_panel(1, 5, k, sar_stretch(imgs{k}), [], ttl{k});
end
save_paper(f, fullfile(OUT, 'fig02_test_images.png'));

fprintf('Saved %s\n', fullfile(OUT, 'fig02_test_images.png'));
fprintf('Image sizes: %s\n', mat2str(size(imgs{1})));
