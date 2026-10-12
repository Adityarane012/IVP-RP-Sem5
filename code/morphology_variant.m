%% ========================================================================
%  morphology_variant.m  —  Task 8: morphological post-processing
%
%  The base paper lists morphological post-processing as optional. Applied
%  to the fused ROSED edge map, after LRT localisation:
%
%    none         no post-processing
%    open-5/10/20 area opening: remove 8-connected components smaller than
%                 5, 10 or 20 px (bwareaopen). Removes isolated false
%                 detections while keeping long edges.
%    close+open10 closing with a 3x3 plus-shaped element, then area opening
%                 at 10 px. Closing bridges one-pixel gaps in LRT edge lines.
%
%  Why area opening and not ordinary opening: after the LRT, boundary edges
%  are one pixel wide, and a standard opening with any element wider than one
%  pixel would erase them entirely.
%
%  Evaluated for fixed-14 and intersection on Fig 6a-d and Fig 6e.
%  Requires fig6abcd.mat, fig6e_images.mat, calibration.mat.
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

LRT_N = 14;  guard = 2;
TAB = '../results/tables';
C = load('calibration.mat');
D = load('../data/simulated/fig6abcd.mat');
E = load('../data/simulated/fig6e_images.mat');
VAR = {'fixed-14', 'intersection'};
MOR = {'none', 'open-5', 'open-10', 'open-20', 'close+open10'};
SE  = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);

hdr = {'Image','Variant','Morph','FN','FP','TPR'};
rows = cell(0, 6);
T = {D.Pa, D.gtA, 'Fig6a'; D.Pb, D.gtB, 'Fig6b'; D.Pc, D.gtC, 'Fig6c'};
fprintf('%-6s %-13s %-13s %5s %6s %6s\n', 'image', 'variant', 'morph', 'FN', 'FP', 'TPR');
for t = 1:3
    S = build_scale_stack(T{t,1}, C.scales, C.taus, C.alpha_pos, C.alpha_neg);
    cf = cfar_fused(S);
    for v = 1:numel(VAR)
        final = fuse_variant(T{t,1}, S, VAR{v}, LRT_N, guard, cf);
        for m = 1:numel(MOR)
            [FN, FP, TPR] = evaluate_perrow(morph(final, MOR{m}, SE), T{t,2}, 2);
            fprintf('%-6s %-13s %-13s %5d %6d %6.1f\n', T{t,3}, VAR{v}, MOR{m}, FN, FP, TPR);
            rows(end+1, :) = {T{t,3}, VAR{v}, MOR{m}, FN, FP, TPR}; %#ok<SAGROW>
        end
    end
end
write_csv(fullfile(TAB, 'morphology_6abc.csv'), hdr, rows);

hdrD = {'Variant','Morph','BoundaryTPR','BrightTPR','DarkTPR','FP'};
rowsD = cell(0, 6);
S = build_scale_stack(D.Pd, C.scales, C.taus, C.alpha_pos, C.alpha_neg);
cf = cfar_fused(S);
fprintf('\n=== Fig6d ===\n');
for v = 1:numel(VAR)
    final = fuse_variant(D.Pd, S, VAR{v}, LRT_N, guard, cf);
    for m = 1:numel(MOR)
        Rm = evaluate_multi(morph(final, MOR{m}, SE), D.gtD, 2);
        fprintf('%-13s %-13s bound %5.1f bright %5.1f dark %5.1f FP %5d\n', VAR{v}, MOR{m}, ...
                Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total);
        rowsD(end+1, :) = {VAR{v}, MOR{m}, Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total}; %#ok<SAGROW>
    end
end
write_csv(fullfile(TAB, 'morphology_6d.csv'), hdrD, rowsD);

hdrE = {'GapPx','Variant','Morph','FP','BothResolvedPct'};
rowsE = cell(0, 5);
fprintf('\n=== Fig6e: both resolved %% / FP ===\n');
[gs, o] = sort(E.GEN_ORDER);
for j = 1:numel(gs)
    i = o(j);  P = E.imgs{i};  gt = E.gts{i};
    S = build_scale_stack(P, C.scales, C.taus, C.alpha_pos, C.alpha_neg);
    cf = cfar_fused(S);
    for v = 1:numel(VAR)
        final = fuse_variant(P, S, VAR{v}, LRT_N, guard, cf);
        fprintf('gap %2d %-13s', gs(j), VAR{v});
        for m = 1:numel(MOR)
            R = evaluate_pair(morph(final, MOR{m}, SE), gt, 2);
            fprintf('  %s %5.1f/%-4d', MOR{m}, R.both, R.FP);
            rowsE(end+1, :) = {gs(j), VAR{v}, MOR{m}, R.FP, R.both}; %#ok<SAGROW>
        end
        fprintf('\n');
    end
end
write_csv(fullfile(TAB, 'morphology_6e.csv'), hdrE, rowsE);
fprintf('\nSaved morphology_6abc.csv, morphology_6d.csv, morphology_6e.csv\n');

function out = morph(E, name, SE)
    switch name
        case 'none',         out = E;
        case 'open-5',       out = bwareaopen(E, 5, 8);
        case 'open-10',      out = bwareaopen(E, 10, 8);
        case 'open-20',      out = bwareaopen(E, 20, 8);
        case 'close+open10', out = bwareaopen(imclose(E, SE), 10, 8);
    end
end
