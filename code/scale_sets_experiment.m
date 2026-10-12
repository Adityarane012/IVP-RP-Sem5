%% ========================================================================
%  scale_sets_experiment.m  —  Which scales should intersection combine?
%
%  Two questions:
%   1. Is {8,14,20} a reasonable choice, or would another set do better?
%   2. Does intersection still resolve close boundaries WITHOUT scale 8?
%      This matters for real data: at P_FA = 1e-3 the s = 8 ROM channel on
%      8-bit land clutter only fires at ~50:1 contrast (PROJECT.md F15).
%
%  Scales 6 and 10 are calibrated here with the same Monte Carlo method as
%  rosed_pipeline.m and saved to calibration_extra.mat. Scales 8/14/20 use
%  calibration.mat unchanged, so every other result stays consistent.
%
%  CFAR is fused over the original scales {8,14,20} for every set, so the
%  only thing that varies between sets is the ROM intersection.
%
%  Uses the same seeds as fig6e_robustness.m (1000*seed + gap), so the
%  images are identical across the two experiments.
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

N_SEEDS = 10;
GAPS    = [6 8 10 12 14 16 24];
SETS    = { [8 14 20], [10 14 20], [14 20], [6 10 14], [8 14], [8 20], [6 14 20] };
SINGLE  = [6 8 10 14 20];
sz = 256;  LRT_N = 14;  guard = 2;  P_FA = 1e-3;
TAB = '../results/tables';  FIG = '../results/figures';

C = load('calibration.mat');

%% Calibrate the extra scales (6, 10)
if exist('calibration_extra.mat', 'file') == 2
    X = load('calibration_extra.mat');
    fprintf('loaded calibration_extra.mat\n');
else
    rng(7);
    X.scales = [6 10];
    X.taus = zeros(1, 2);  X.alpha_pos = zeros(1, 2);  X.alpha_neg = zeros(1, 2);
    for i = 1:2
        s = X.scales(i);
        X.taus(i) = calibrate_rom(s, P_FA, 3, sz);
        N = s^2;  k = round(0.75*N);
        X.alpha_pos(i) = solve_alpha(N, k, P_FA);
        X.alpha_neg(i) = solve_alpha(N, k, 1 - P_FA);
        [rp, rn] = verify_cfar(s, X.alpha_pos(i), X.alpha_neg(i), sz);
        fprintf('scale %2d: tau %.4f  alpha %.3f  gate pos %.2e neg %.2e\n', ...
                s, X.taus(i), X.alpha_pos(i), rp, rn);
        if rp < P_FA/10 || rp > P_FA*10
            error('CFAR gate failed for extra scale %d', s);
        end
    end
    save('calibration_extra.mat', '-struct', 'X', '-v7');
    save(fullfile(TAB, 'calibration_extra.mat'), '-struct', 'X', '-v7');
end

% combined calibration over all scales used
allS = [C.scales, X.scales];
allT = [C.taus, X.taus];
allP = [C.alpha_pos, X.alpha_pos];
allN = [C.alpha_neg, X.alpha_neg];
[allS, o] = sort(allS);  allT = allT(o);  allP = allP(o);  allN = allN(o);
fprintf('scales %s  tau %s\n', mat2str(allS), mat2str(allT, 4));

nSet = numel(SETS);  nSgl = numel(SINGLE);  nG = numel(GAPS);
nCol = nSet + nSgl;
% res(seed, gap, column, metric), metric = [both, FP]
res = zeros(N_SEEDS, nG, nCol, 2);
orig = ismember(allS, [8 14 20]);

t0 = tic;
parfor s = 1:N_SEEDS
    r = zeros(nG, nCol, 2);
    for gi = 1:nG
        g = GAPS(gi);
        rng(1000*s + g);
        [P, gt] = make_fig6e(sz, g);
        Sall = build_scale_stack(P, allS, allT, allP, allN);
        cf   = cfar_fused(Sall(orig));
        for c = 1:nSet
            idx = ismember(allS, SETS{c});
            E = select_intersection(Sall(idx));
            final = localize_rom(P, E, LRT_N, guard) | cf;
            R = evaluate_pair(final, gt, 2);
            r(gi, c, :) = [R.both, R.FP];
        end
        for c = 1:nSgl
            E = Sall(allS == SINGLE(c)).rom;
            final = localize_rom(P, E, LRT_N, guard) | cf;
            R = evaluate_pair(final, gt, 2);
            r(gi, nSet + c, :) = [R.both, R.FP];
        end
    end
    res(s, :, :, :) = r;
end
fprintf('done in %.0f s\n', toc(t0));

%% Table
names = [cellfun(@(x) ['AND' mat2str(x)], SETS, 'UniformOutput', false), ...
         arrayfun(@(x) sprintf('fixed-%d', x), SINGLE, 'UniformOutput', false)];
hdr = {'Set','GapPx','BothMean','BothStd','FPMean','FPStd'};
rows = cell(0, 6);
fprintf('\n%-16s', 'set \\ gap');  fprintf('%9d', GAPS);  fprintf('%11s\n', 'mean FP');
for c = 1:nCol
    fprintf('%-16s', names{c});
    for gi = 1:nG
        b = res(:, gi, c, 1);  fp = res(:, gi, c, 2);
        fprintf('%9.1f', mean(b));
        rows(end+1, :) = {names{c}, GAPS(gi), mean(b), std(b), mean(fp), std(fp)}; %#ok<SAGROW>
    end
    fprintf('%11.1f\n', mean(reshape(res(:, :, c, 2), [], 1)));
end
fprintf('(cells: both resolved %%, mean over %d seeds)\n', N_SEEDS);
write_csv(fullfile(TAB, 'scale_sets_summary.csv'), hdr, rows);

%% Summary scores
fprintf('\n%-16s %14s %14s %12s\n', 'set', 'mean both%', 'min both%', 'mean FP');
for c = 1:nCol
    b = res(:, :, c, 1);  fp = res(:, :, c, 2);
    fprintf('%-16s %14.1f %14.1f %12.1f\n', names{c}, mean(b(:)), min(mean(b, 1)), mean(fp(:)));
end
fprintf('\nmin both%% = worst gap after averaging over seeds.\n');
fprintf('Sets without 8 px ([10 14 20], [14 20]) answer the real-data question.\n');

%% Figure
f = new_figure('Scale sets');
subplot(1,2,1); hold on;
for c = 1:nSet
    plot(GAPS, mean(res(:, :, c, 1), 1), 'o-', 'LineWidth', 1.2);
end
xlabel('boundary separation (px)'); ylabel('both resolved (% rows)');
title('Intersection sets'); legend(names(1:nSet), 'Location', 'southeast'); grid on;
subplot(1,2,2); hold on;
for c = 1:nSet
    plot(GAPS, mean(res(:, :, c, 2), 1), 'o-', 'LineWidth', 1.2);
end
xlabel('boundary separation (px)'); ylabel('false positives');
title('False positives'); grid on;
save_png(f, fullfile(FIG, 'scale_sets.png'));
fprintf('Saved scale_sets_summary.csv, scale_sets.png\n');
