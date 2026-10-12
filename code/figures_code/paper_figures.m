%% ========================================================================
%  paper_figures.m  —  Regenerate Figs. 8, 9 and 10 of the paper in MATLAB
%
%  WHY THIS EXISTS
%   Fig. 8  was drawn outside MATLAB; it must come from the real pipeline.
%   Fig. 9  was two screenshots with different backgrounds and labels.
%   Fig. 10 had faint grey titles: MATLAB's dark theme keeps text light
%           even when the figure background is set to white.
%
%  WHAT IT USES (no new random numbers, nothing is recomputed differently)
%   calibration.mat                         written by rosed_pipeline.m
%   ../data/simulated/fig6abcd.mat          written by rosed_pipeline.m
%   ../data/simulated/fig6e_images.mat      written by fig6e_extended.m
%   ../results/tables/calibration_real_pooled.mat   real thresholds
%   ../data/crops/*__real_pair*.tif         the real pair crop
%
%  OUTPUT  ../results/figures/paper/fig08_imageE_10px.png
%          ../results/figures/paper/fig09_corroborated_selection.png
%          ../results/figures/paper/fig10_real_pair.png
%
%  CHECK   The numbers printed below must match the paper:
%           Fig. 8  -> Table IX, image E row (seed-42 realisation)
%           Fig. 10 -> Table X, pair row, and the 2.66 / 20.18 / 30.61 %
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

LRT_N = 14;  guard = 2;
SIM = '../data/simulated';  TAB = '../results/tables';
OUT = '../results/figures/paper';
if ~exist(OUT, 'dir'), mkdir(OUT); end
C = load('calibration.mat');
X = char(215);                                   % multiplication sign
if exist('OCTAVE_VERSION', 'builtin'), X = "\xC3\x97"; end   % UTF-8 in Octave

%% ---------------- Fig. 8: image E, 10 px, fixed-14 vs intersection -----
E = load(fullfile(SIM, 'fig6e_images.mat'));
i10 = find(E.GEN_ORDER == 10, 1);
P  = E.imgs{i10};  gt = E.gts{i10};
S  = build_scale_stack(P, C.scales, C.taus, C.alpha_pos, C.alpha_neg);
cf = cfar_fused(S);
Ef = fuse_variant(P, S, 'fixed-14',     LRT_N, guard, cf);
Ei = fuse_variant(P, S, 'intersection', LRT_N, guard, cf);
Rf = evaluate_pair(Ef, gt, 2);
Ri = evaluate_pair(Ei, gt, 2);
fprintf('Fig. 8  fixed-14:     both resolved %.1f%%, FP %d\n', Rf.both, Rf.FP);
fprintf('Fig. 8  intersection: both resolved %.1f%%, FP %d\n', Ri.both, Ri.FP);
fprintf('        (Table IX expects 23.4%% / 228 and 99.6%% / 164)\n');

cols = max(1, gt.edge1(1) - 40) : min(size(P, 2), gt.edge2(2) + 40);
G = sar_stretch(P);
f = paper_fig('Fig 8', [100 100 1100 560]);
show_panel(1, 3, 1, G(:, cols), [],          '(a) Input');
show_panel(1, 3, 2, G(:, cols), Ef(:, cols), sprintf('(b) Fixed 14 %s 14: %.1f%%', X, Rf.both));
show_panel(1, 3, 3, G(:, cols), Ei(:, cols), sprintf('(c) Intersection: %.1f%%', Ri.both));
save_paper(f, fullfile(OUT, 'fig08_imageE_10px.png'));

%% ---------------- Fig. 9: corroborated scale selection ------------------
D  = load(fullfile(SIM, 'fig6abcd.mat'));
SB = build_scale_stack(D.Pb, C.scales, C.taus, C.alpha_pos, C.alpha_neg);
[~, ~, selB] = fuse_variant(D.Pb, SB, 'corroborated', LRT_N, guard, cfar_fused(SB));
[~, ~, selE] = fuse_variant(P,    S,  'corroborated', LRT_N, guard, cf);
fprintf('Fig. 9  image E: share of selected pixels at 20 %s 20 = %.1f%%\n', X, ...
        100 * nnz(selE == 3) / max(nnz(selE), 1));

CMAP = [1.00 1.00 1.00;      % 0 none
        0.85 0.33 0.10;      % 1 = 8 x 8
        0.00 0.45 0.74;      % 2 = 14 x 14
        0.20 0.20 0.20];     % 3 = 20 x 20
f = paper_fig('Fig 9', [100 100 1200 520]);
sel_panel(1, selB, CMAP, '(a) Image B (bright line)');
sel_panel(2, selE, CMAP, '(b) Image E, 10 px separation');
save_paper(f, fullfile(OUT, 'fig09_corroborated_selection.png'));

%% ---------------- Fig. 10: real pair crop, ROM + LRT channel ------------
R  = load(fullfile(TAB, 'calibration_real_pooled.mat'));
fl = dir(fullfile('../data/crops', '*__real_pair*.tif'));
if isempty(fl), error('Real pair crop not found in ../data/crops'); end
A  = double(imread(fullfile(fl(1).folder, fl(1).name)));
SR = build_scale_stack(A, R.SCALES, R.taus_real, R.alphas_real, R.alphas_neg_real);
nocfar = false(size(A));
Rfx = fuse_variant(A.^2, SR, 'fixed-14',     LRT_N, guard, nocfar);
Rin = fuse_variant(A.^2, SR, 'intersection', LRT_N, guard, nocfar);
fprintf('Fig. 10 ROM + LRT flagged: fixed-14 %.2f%%, intersection %.2f%%\n', ...
        100 * mean(Rfx(:)), 100 * mean(Rin(:)));
fprintf('        (Table X expects 3.81%% and 0.69%%)\n');

GR = sar_stretch(A, 'linear');
f = paper_fig('Fig 10', [50 50 1500 560]);
show_panel(1, 3, 1, GR, [],  '(a) Input');
show_panel(1, 3, 2, GR, Rfx, sprintf('(b) ROM + LRT, fixed 14 %s 14', X));
show_panel(1, 3, 3, GR, Rin, '(c) ROM + LRT, intersection');
for i = 1:numel(SR)
    pct = 100 * mean(SR(i).rom(:));
    fprintf('Fig. 10 ROM fires at s = %2d: %.2f%%\n', SR(i).scale, pct);
end
fprintf('        (the text expects 2.66 / 20.18 / 30.61%%)\n');
save_paper(f, fullfile(OUT, 'fig10_real_pair.png'));

fprintf('\nSaved three figures to %s\n', OUT);
% Helpers: paper_fig.m, show_panel.m, sel_panel.m, style_title.m, save_paper.m
