%% ========================================================================
%  fig6e_robustness.m  —  Is the Fig6e result robust to the noise draw?
%
%  fig6e_extended.m uses one speckle realisation per gap. This repeats the
%  experiment over N_SEEDS independent realisations and reports mean and
%  standard deviation, so the paper's headline numbers are not one draw.
%
%  Seed for (seed s, gap g) is 1000*s + g: deterministic, independent of
%  run order, and safe under parfor.
%
%  Also tests the claim directly: does intersection give fewer false
%  positives than EVERY single scale while matching the best resolution?
%
%  Requires calibration.mat in code/. Uses parfor if Parallel Computing
%  Toolbox is available; otherwise runs serially with identical results.
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

N_SEEDS = 10;
GAPS    = [6 8 10 12 14 16 24];
VAR     = {'fixed-14','union','intersection','corroborated','fixed-8','fixed-20'};
sz = 256;  LRT_N = 14;  guard = 2;
TAB = '../results/tables';  FIG = '../results/figures';

C = load('calibration.mat');
scales = C.scales;  taus = C.taus;  ap = C.alpha_pos;  an = C.alpha_neg;
nG = numel(GAPS);  nV = numel(VAR);

% res(seed, gap, variant, metric), metric = [both, FP, e1, e2, twoplus]
res = zeros(N_SEEDS, nG, nV, 5);
t0 = tic;
parfor s = 1:N_SEEDS
    r = zeros(nG, nV, 5);
    for gi = 1:nG
        g = GAPS(gi);
        rng(1000*s + g);
        [P, gt] = make_fig6e(sz, g);
        S  = build_scale_stack(P, scales, taus, ap, an);
        cf = cfar_fused(S);
        lo = max(1, min(gt.edge1) - 20);  hi = min(sz, max(gt.edge2) + 20);
        for v = 1:nV
            [final, romsel] = fuse_variant(P, S, VAR{v}, LRT_N, guard, cf);
            R = evaluate_pair(final, gt, 2);
            W = run_stats(romsel, lo, hi);
            r(gi, v, :) = [R.both, R.FP, R.e1, R.e2, W.frac_two_plus];
        end
    end
    res(s, :, :, :) = r;
end
fprintf('%d seeds x %d gaps x %d variants in %.0f s\n', N_SEEDS, nG, nV, toc(t0));

%% Summary table
hdr = {'GapPx','Variant','BothMean','BothStd','FPMean','FPStd','Edge1Mean','Edge2Mean','TwoPlusRunsMean','NSeeds'};
rows = cell(0, numel(hdr));
fprintf('\n%-6s %-14s %16s %16s\n', 'gap', 'variant', 'both%% mean±std', 'FP mean±std');
for gi = 1:nG
    for v = 1:nV
        b  = res(:, gi, v, 1);  fp = res(:, gi, v, 2);
        fprintf('%-6d %-14s %8.1f ± %5.1f %8.1f ± %5.1f\n', GAPS(gi), VAR{v}, ...
                mean(b), std(b), mean(fp), std(fp));
        rows(end+1, :) = {GAPS(gi), VAR{v}, mean(b), std(b), mean(fp), std(fp), ...
                          mean(res(:,gi,v,3)), mean(res(:,gi,v,4)), mean(res(:,gi,v,5)), N_SEEDS}; %#ok<SAGROW>
    end
    fprintf('\n');
end
write_csv(fullfile(TAB, 'fig6e_robustness_summary.csv'), hdr, rows);

% raw per-seed results, so any number can be re-derived
hraw = {'Seed','GapPx','Variant','Both','FP','Edge1','Edge2','TwoPlusRuns'};
raw = cell(0, numel(hraw));
for s = 1:N_SEEDS
    for gi = 1:nG
        for v = 1:nV
            m = squeeze(res(s, gi, v, :))';
            raw(end+1, :) = {s, GAPS(gi), VAR{v}, m(1), m(2), m(3), m(4), m(5)}; %#ok<SAGROW>
        end
    end
end
write_csv(fullfile(TAB, 'fig6e_robustness_raw.csv'), hraw, raw);

%% Claim tests, counted over every (seed, gap) pair
iI = find(strcmp(VAR, 'intersection'));
iF = find(strcmp(VAR, 'fixed-14'));
single = find(strncmp(VAR, 'fixed-', 6));
bothI = res(:, :, iI, 1);  fpI = res(:, :, iI, 2);
bestSingleBoth = max(res(:, :, single, 1), [], 3);
minSingleFP    = min(res(:, :, single, 2), [], 3);
n = numel(bothI);
fprintf('=== Claim tests over %d (seed, gap) pairs ===\n', n);
fprintf('intersection FP < every single scale          : %3d / %d\n', nnz(fpI < minSingleFP), n);
fprintf('intersection both%% >= best single scale - 1   : %3d / %d\n', nnz(bothI >= bestSingleBoth - 1), n);
fprintf('intersection both%% > fixed-14 + 1             : %3d / %d\n', nnz(bothI > res(:,:,iF,1) + 1), n);
fprintf('intersection both%% < fixed-14 - 1             : %3d / %d\n', nnz(bothI < res(:,:,iF,1) - 1), n);
tp = res(:, :, :, 5);  bb = res(:, :, :, 1);
cc = corrcoef(tp(:), bb(:));
fprintf('corr(rows with 2+ runs, both resolved), all cases: r = %.3f (n = %d)\n', cc(1,2), numel(tp));

%% Figure: mean ± std of both-resolved vs gap
f = new_figure('Fig6e robustness');  hold on;
mk = {'o-','s-','d-','^-','v--','x--'};
for v = 1:nV
    mu = mean(res(:, :, v, 1), 1);  sd = std(res(:, :, v, 1), 0, 1);
    h = errorbar(GAPS, mu, sd, mk{v});
    set(h, 'LineWidth', 1.2);
end
xlabel('boundary separation (px)'); ylabel('both resolved (% rows), mean ± std');
title(sprintf('Fig6e over %d speckle realisations', N_SEEDS));
legend(VAR, 'Location', 'southeast'); grid on;
save_png(f, fullfile(FIG, 'fig6e_robustness.png'));
fprintf('\nSaved fig6e_robustness_summary.csv, fig6e_robustness_raw.csv, fig6e_robustness.png\n');
