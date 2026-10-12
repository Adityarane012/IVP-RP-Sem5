%% ========================================================================
%  classical_baselines.m  —  Task 4 and LoG/Canny baselines
%
%  Uses MATLAB's edge() with its default automatic thresholds:
%    sobel  (thinned, the standard implementation)
%    log    (Laplacian of Gaussian)
%    canny
%  plus our unthinned Sobel (sobel_edge.m), which produced the 3984 figure
%  in the Introduction, so both can be reported side by side.
%
%  Input is the normalised amplitude, as in rosed_pipeline.m Section 6.
%  Same evaluation code as ROSED, so the numbers are directly comparable.
%
%  Requires ../data/simulated/fig6abcd.mat and fig6e_images.mat.
%  Needs Image Processing Toolbox (edge).
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

TAB = '../results/tables';
D = load('../data/simulated/fig6abcd.mat');
E = load('../data/simulated/fig6e_images.mat');
M = {'sobel-unthinned', 'sobel', 'log', 'canny'};

hdr  = {'Method','Image','FN','FP','TPR'};
rows = cell(0, 5);
fprintf('=== Fig 6a/b/c ===\n%-16s %-6s %6s %8s %7s\n', 'method', 'image', 'FN', 'FP', 'TPR%');
T = {D.Pa, D.gtA, 'Fig6a'; D.Pb, D.gtB, 'Fig6b'; D.Pc, D.gtC, 'Fig6c'};
for t = 1:3
    for m = 1:numel(M)
        e = run_method(T{t,1}, M{m});
        [FN, FP, TPR] = evaluate_perrow(e, T{t,2}, 2);
        fprintf('%-16s %-6s %6d %8d %7.1f\n', M{m}, T{t,3}, FN, FP, TPR);
        rows(end+1, :) = {M{m}, T{t,3}, FN, FP, TPR}; %#ok<SAGROW>
    end
end
write_csv(fullfile(TAB, 'classical_6abc.csv'), hdr, rows);

hdrD = {'Method','BoundaryTPR','BrightTPR','DarkTPR','FP'};
rowsD = cell(0, 5);
fprintf('\n=== Fig 6d ===\n');
for m = 1:numel(M)
    Rm = evaluate_multi(run_method(D.Pd, M{m}), D.gtD, 2);
    fprintf('%-16s bound %5.1f  bright %5.1f  dark %5.1f  FP %6d\n', M{m}, ...
            Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total);
    rowsD(end+1, :) = {M{m}, Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total}; %#ok<SAGROW>
end
write_csv(fullfile(TAB, 'classical_6d.csv'), hdrD, rowsD);

hdrE = {'Method','GapPx','Edge1TPR','Edge2TPR','FP','BothResolvedPct'};
rowsE = cell(0, 6);
fprintf('\n=== Fig 6e (single realisation): both resolved %% / FP ===\n%-16s', 'method \\ gap');
[gs, o] = sort(E.GEN_ORDER);
fprintf('%13d', gs); fprintf('\n');
for m = 1:numel(M)
    fprintf('%-16s', M{m});
    for j = 1:numel(gs)
        i = o(j);
        R = evaluate_pair(run_method(E.imgs{i}, M{m}), E.gts{i}, 2);
        fprintf('%7.1f/%-5d', R.both, R.FP);
        rowsE(end+1, :) = {M{m}, gs(j), R.e1, R.e2, R.FP, R.both}; %#ok<SAGROW>
    end
    fprintf('\n');
end
write_csv(fullfile(TAB, 'classical_6e.csv'), hdrE, rowsE);
fprintf('\nSaved classical_6abc.csv, classical_6d.csv, classical_6e.csv\n');
fprintf('\nREAD TPR ONLY TOGETHER WITH FP. A method that marks most of the image\n');
fprintf('scores high TPR and high "both resolved" by flooding, not by detecting.\n');
fprintf('Update the Introduction: use the thinned sobel FP for Fig6a.\n');

function e = run_method(P, m)
    import rosed.*          % scoped per function, so repeated here
    I = norm01(sqrt(P));
    switch m
        case 'sobel-unthinned', e = sobel_edge(I);
        case 'sobel',           e = edge(I, 'sobel');
        case 'log',             e = edge(I, 'log');
        case 'canny',           e = edge(I, 'canny');
    end
end
