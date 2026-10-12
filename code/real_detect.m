%% ========================================================================
%  real_detect.m  —  Task 3: ROSED on real UMBRA crops (faculty demo)
%
%  Runs every variant on each ../data/crops/*__real_*.tif with the
%  empirically calibrated real-data thresholds, and saves:
%    - a three-panel figure per crop: input | fixed-14 | intersection
%    - a per-scale ROM figure, showing how often each scale fires
%    - real_detect_summary.csv
%
%  DOMAINS (see PROJECT.md F2, F14)
%    ROM and OS-CFAR on AMPLITUDE: the uint8 values as supplied, which is
%    what the thresholds were calibrated on.
%    LRT localisation on INTENSITY (A.^2): the LRT picks the most likely
%    split position under an exponential model; it applies no threshold.
%
%  There is no ground truth for real scenes, so the summary reports what can
%  be measured without it: fraction of pixels flagged, number of connected
%  edge segments, and mean segment length. Longer, fewer segments indicate
%  coherent edges; many tiny segments indicate scattered false alarms.
%
%  Requires calibration_real_pooled.mat with alphas_neg_real
%  (run real_calibrate_neg.m first).
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

LRT_N = 14;  guard = 2;
TAB = '../results/tables';  FIG = '../results/figures/real';
R = load(fullfile(TAB, 'calibration_real_pooled.mat'));
if ~isfield(R, 'alphas_neg_real')
    error('alphas_neg_real missing. Run real_calibrate_neg.m first.');
end
SCALES = R.SCALES;
fprintf('real thresholds  tau %s  alpha_pos %s  alpha_neg %s\n', ...
        mat2str(R.taus_real, 4), mat2str(R.alphas_real, 4), mat2str(R.alphas_neg_real, 4));

files = dir(fullfile('../data/crops', '*__real_*.tif'));
if isempty(files)
    error('No *__real_*.tif crops found. Add real_ rows to crops.csv and extract.');
end

VAR = {'fixed-14', 'intersection', 'union', 'corroborated', 'fixed-8', 'fixed-20'};
SE  = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);
hdr = {'Crop','Variant','Morph','PctFlagged','Segments','MeanSegmentPx','RomPct8','RomPct14','RomPct20'};
rows = cell(0, numel(hdr));

for fi = 1:numel(files)
    [~, name] = fileparts(files(fi).name);
    short = regexprep(name, '^.*__', '');
    A = double(imread(fullfile(files(fi).folder, files(fi).name)));
    fprintf('\n=== %s  (%dx%d, mean %.1f, CV %.3f) ===\n', short, size(A,1), size(A,2), ...
            mean(A(:)), std(A(:))/mean(A(:)));

    S  = build_scale_stack(A, SCALES, R.taus_real, R.alphas_real, R.alphas_neg_real);
    cf = cfar_fused(S);
    romPct = arrayfun(@(x) 100*mean(x.rom(:)), S);
    fprintf('ROM fires on %.2f%% / %.2f%% / %.2f%% of pixels at s = %s\n', romPct, mat2str(SCALES));
    fprintf('CFAR (fused) flags %.2f%% of pixels\n', 100*mean(cf(:)));

    fprintf('%-13s %-13s %8s %9s %10s\n', 'variant', 'morph', 'flag%', 'segments', 'mean px');
    maps = struct();
    for v = 1:numel(VAR)
        final = fuse_variant(A.^2, S, VAR{v}, LRT_N, guard, cf);
        for mm = 1:2
            if mm == 1
                E = final;  mname = 'none';
            else
                E = bwareaopen(imclose(final, SE), 10, 8);  mname = 'close+open10';
            end
            [nseg, meanlen] = segment_stats(E);
            fprintf('%-13s %-13s %8.2f %9d %10.1f\n', VAR{v}, mname, 100*mean(E(:)), nseg, meanlen);
            rows(end+1, :) = {short, VAR{v}, mname, 100*mean(E(:)), nseg, meanlen, ...
                              romPct(1), romPct(2), romPct(3)}; %#ok<SAGROW>
            if mm == 1 && strcmp(VAR{v}, 'fixed-14'),     maps.fixed = E; end
            if mm == 1 && strcmp(VAR{v}, 'intersection'), maps.inter = E; end
        end
    end

    % three-panel figure
    G = sar_stretch(A, 'linear');
    f = new_figure(['Real — ' short]);
    set(f, 'Position', [50 50 1500 540]);
    subplot(1,3,1); imshow(G); title(sprintf('(a) %s, input', strrep(short, '_', ' ')));
    subplot(1,3,2); imshow(overlay(G, maps.fixed)); title('(b) fixed 14 x 14');
    subplot(1,3,3); imshow(overlay(G, maps.inter)); title('(c) intersection {8,14,20}');
    save_png(f, fullfile(FIG, [short '_3panel.png']));

    % per-scale ROM figure: shows the s = 8 gating directly
    f = new_figure(['Real ROM per scale — ' short]);
    set(f, 'Position', [50 50 1500 540]);
    for i = 1:numel(S)
        subplot(1, numel(S), i); imshow(overlay(G, S(i).rom));
        title(sprintf('ROM s = %d  (%.2f%%)', S(i).scale, romPct(i)));
    end
    save_png(f, fullfile(FIG, [short '_rom_per_scale.png']));
end

write_csv(fullfile(TAB, 'real_detect_summary.csv'), hdr, rows);
fprintf('\nSaved real_detect_summary.csv and figures in %s\n', FIG);
fprintf('Expect intersection to be sparse: at s = 8 the real threshold only\n');
fprintf('passes ~50:1 contrast (PROJECT.md F15). Report it as observed.\n');

function [n, meanlen] = segment_stats(E)
    L = bwlabel(E, 8);
    n = max(L(:));
    if n == 0
        meanlen = 0;
    else
        meanlen = nnz(E) / n;
    end
end

function rgb = overlay(G, E)
    r = G;  g = G;  b = G;
    r(E) = 1;  g(E) = 0;  b(E) = 0;
    rgb = cat(3, r, g, b);
end
