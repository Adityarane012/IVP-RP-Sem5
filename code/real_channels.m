%% ========================================================================
%  real_channels.m  —  Separate the ROM and CFAR channels on real crops
%
%  WHY THIS EXISTS
%  In real_detect.m the fused maps for fixed-14 and intersection look almost
%  identical. The CFAR channel is the same for every variant and it
%  dominates, so it hides the ROM + LRT channel, which is the only place the
%  variants differ. This script shows that channel on its own.
%
%  It also measures a real-data finding: on 8-bit UMBRA products, a large
%  share of pixels sit at exactly 1 (the provider's floor). The CFAR false-
%  alarm rate then depends strongly on local brightness, so it is no longer
%  constant. The table reports it per brightness band.
%
%  Outputs
%    ..\results\figures\real\<crop>_rom_channel.png
%    ..\results\tables\real_cfar_by_brightness.csv
%    ..\results\tables\real_rom_channel.csv
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

LRT_N = 14;  guard = 2;  S14 = 2;             % index of scale 14
TAB = '../results/tables';  FIG = '../results/figures/real';
R = load(fullfile(TAB, 'calibration_real_pooled.mat'));
if ~isfield(R, 'alphas_neg_real'), error('Run real_calibrate_neg.m first.'); end
SCALES = R.SCALES;
files = dir(fullfile('../data/crops', '*__real_*.tif'));
if isempty(files), error('No *__real_*.tif crops found.'); end

BANDS = [0 15 25 35 50 75 256];                % local-mean bands (uint8 units)
hdrB = {'Crop','BandLo','BandHi','SharePct','PosCFARPct','NegCFARPct'};
rowsB = cell(0, 6);
hdrR = {'Crop','Variant','Morph','PctFlagged','Segments','MeanSegmentPx'};
rowsR = cell(0, 6);
SE = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);

for fi = 1:numel(files)
    [~, name] = fileparts(files(fi).name);
    short = regexprep(name, '^.*__', '');
    A = double(imread(fullfile(files(fi).folder, files(fi).name)));
    fprintf('\n=== %s ===\n', short);
    fprintf('pixels exactly at the floor value %d: %.1f%%\n', min(A(:)), 100*mean(A(:) == min(A(:))));

    S = build_scale_stack(A, SCALES, R.taus_real, R.alphas_real, R.alphas_neg_real);

    % ---- CFAR rate by local brightness, s = 14 ---------------------------
    s = SCALES(S14);  b = s;
    loc = conv2(A, ones(s) / s^2, 'same');
    loc = loc(b:end-b, b:end-b);
    cp  = S(S14).cfar_p(b:end-b, b:end-b);
    cn  = S(S14).cfar_n(b:end-b, b:end-b);
    fprintf('CFAR at s = 14: positive %.2f%%, negative %.2f%% (target 0.10%%)\n', ...
            100*mean(cp(:)), 100*mean(cn(:)));
    fprintf('%12s %10s %10s %10s\n', 'local mean', 'share', 'pos CFAR', 'neg CFAR');
    for j = 1:numel(BANDS) - 1
        m = loc >= BANDS(j) & loc < BANDS(j+1);
        if ~any(m(:)), continue; end
        fprintf('%5d-%-6d %9.1f%% %9.2f%% %9.2f%%\n', BANDS(j), BANDS(j+1), ...
                100*mean(m(:)), 100*mean(cp(m)), 100*mean(cn(m)));
        rowsB(end+1, :) = {short, BANDS(j), BANDS(j+1), 100*mean(m(:)), ...
                           100*mean(cp(m)), 100*mean(cn(m))}; %#ok<SAGROW>
    end

    % ---- ROM + LRT channel only, no CFAR ---------------------------------
    nocfar = false(size(A));
    VAR = {'fixed-14', 'intersection', 'fixed-8', 'fixed-20'};
    maps = struct();
    fprintf('ROM + LRT channel only:\n%-13s %-13s %8s %9s %10s\n', ...
            'variant', 'morph', 'flag%', 'segments', 'mean px');
    for v = 1:numel(VAR)
        E0 = fuse_variant(A.^2, S, VAR{v}, LRT_N, guard, nocfar);
        for mm = 1:2
            if mm == 1, E = E0; mname = 'none';
            else,       E = bwareaopen(imclose(E0, SE), 10, 8); mname = 'close+open10'; end
            L = bwlabel(E, 8);  n = max(L(:));
            ml = 0; if n > 0, ml = nnz(E) / n; end
            fprintf('%-13s %-13s %8.2f %9d %10.1f\n', VAR{v}, mname, 100*mean(E(:)), n, ml);
            rowsR(end+1, :) = {short, VAR{v}, mname, 100*mean(E(:)), n, ml}; %#ok<SAGROW>
        end
        if strcmp(VAR{v}, 'fixed-14'),     maps.fixed = E0; end
        if strcmp(VAR{v}, 'intersection'), maps.inter = E0; end
    end

    G = sar_stretch(A, 'linear');
    f = new_figure(['ROM channel — ' short]);
    set(f, 'Position', [50 50 1500 540]);
    subplot(1,3,1); imshow(G); title(sprintf('(a) %s, input', strrep(short, '_', ' ')));
    subplot(1,3,2); imshow(overlay(G, maps.fixed)); title('(b) ROM + LRT, fixed 14 x 14');
    subplot(1,3,3); imshow(overlay(G, maps.inter)); title('(c) ROM + LRT, intersection');
    save_png(f, fullfile(FIG, [short '_rom_channel.png']));
end

write_csv(fullfile(TAB, 'real_cfar_by_brightness.csv'), hdrB, rowsB);
write_csv(fullfile(TAB, 'real_rom_channel.csv'), hdrR, rowsR);
fprintf('\nSaved real_cfar_by_brightness.csv, real_rom_channel.csv and *_rom_channel.png\n');

function rgb = overlay(G, E)
    r = G; g = G; b = G;
    r(E) = 1; g(E) = 0; b(E) = 0;
    rgb = cat(3, r, g, b);
end
