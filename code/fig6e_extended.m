%% ========================================================================
%  fig6e_extended.m  —  Tasks 5 and 6
%
%  5. Finer gap sweep: 6, 8, 10, 12, 14, 16, 24 px
%  6. Run-width diagnostic: measures the width and number of detection runs
%     reaching the LRT, to evidence the mechanism directly.
%
%  Images are generated in the order 6, 10, 16, 24, 8, 12, 14 after
%  rng(42). The first four are therefore identical to fig6e_experiment.m,
%  and their results are regression-checked against results_fig6e.csv.
%
%  Also adds single-scale variants fixed-8 and fixed-20, so the table shows
%  what each scale does on its own next to the combinations.
%
%  Requires calibration.mat in code/ (written by rosed_pipeline.m).
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m
rng(42);

sz = 256;  LRT_N = 14;  guard = 2;
TAB = '../results/tables';  SIM = '../data/simulated';  FIG = '../results/figures';

C = load('calibration.mat');
scales = C.scales;
fprintf('calibration: tau = %s\n', mat2str(C.taus, 4));

GEN_ORDER = [6 10 16 24 8 12 14];     % do not reorder: preserves v1 images
GAPS      = sort(GEN_ORDER);
imgs = cell(1, numel(GEN_ORDER));  gts = cell(1, numel(GEN_ORDER));
for i = 1:numel(GEN_ORDER)
    [imgs{i}, gts{i}] = make_fig6e(sz, GEN_ORDER(i));
end
save(fullfile(SIM, 'fig6e_images.mat'), 'GEN_ORDER', 'imgs', 'gts', 'sz', '-v7');

VAR = {'fixed-14', 'union', 'intersection', 'corroborated', 'fixed-8', 'fixed-20'};
hdrR = {'GapPx', 'Variant', 'Edge1TPR', 'Edge2TPR', 'FP', 'BothResolvedPct'};
hdrW = {'GapPx', 'Variant', 'MeanRunWidth', 'MeanRunsPerRow', 'RowsTwoPlusRunsPct', 'BothResolvedPct'};
rowsR = cell(0, 6);  rowsW = cell(0, 6);

for g = GAPS
    i = find(GEN_ORDER == g);
    P = imgs{i};  gt = gts{i};
    S  = build_scale_stack(P, scales, C.taus, C.alpha_pos, C.alpha_neg);
    cf = cfar_fused(S);

    % analysis window: 20 px either side of the boundary pair
    lo = max(1, min(gt.edge1) - 20);
    hi = min(sz, max(gt.edge2) + 20);

    fprintf('\n--- gap = %d px ---\n', g);
    fprintf('%-14s %9s %9s %6s %9s | %8s %8s %8s\n', 'variant', 'edge1', 'edge2', ...
            'FP', 'both%', 'runW', 'runs/row', '2+runs%');
    for v = 1:numel(VAR)
        [final, romsel] = fuse_variant(P, S, VAR{v}, LRT_N, guard, cf);
        R  = evaluate_pair(final, gt, 2);
        W  = run_stats(romsel, lo, hi);
        fprintf('%-14s %9.1f %9.1f %6d %9.1f | %8.1f %8.2f %8.1f\n', VAR{v}, ...
                R.e1, R.e2, R.FP, R.both, W.mean_width, W.mean_runs, W.frac_two_plus);
        rowsR(end+1, :) = {g, VAR{v}, R.e1, R.e2, R.FP, R.both}; %#ok<SAGROW>
        rowsW(end+1, :) = {g, VAR{v}, W.mean_width, W.mean_runs, W.frac_two_plus, R.both}; %#ok<SAGROW>
    end
end

%% Regression: original four gaps, original four variants
fprintf('\n=== Regression against fig6e_experiment.m (gaps 6/10/16/24) ===\n');
orig = ismember(cell2mat(rowsR(:,1)), [6 10 16 24]) & ...
       ismember(rowsR(:,2), {'fixed-14','union','intersection','corroborated'});
ok = compare_to_baseline({fullfile(TAB,'results_fig6e.csv'), 'results_fig6e.csv'}, ...
                         hdrR, rowsR(orig, :), {'GapPx','Variant'}, ...
                         {'Edge1TPR','Edge2TPR','FP','BothResolvedPct'}, 1e-3);
REGRESSION_OK = ok;
if ~ok
    fprintf('REGRESSION FAILED: the original gaps no longer reproduce. Stop here.\n');
end

write_csv(fullfile(TAB, 'results_fig6e_extended.csv'), hdrR, rowsR);
write_csv(fullfile(TAB, 'run_width_fig6e.csv'), hdrW, rowsW);

%% Does run structure predict resolution?
two  = cell2mat(rowsW(:,5));  both = cell2mat(rowsW(:,6));
cc = corrcoef(two, both);
fprintf('\nCorrelation, rows with 2+ runs vs both resolved: r = %.3f (n = %d)\n', ...
        cc(1,2), numel(two));
fprintf('A high r means resolution is explained by run separation before the LRT.\n');

%% Figures
f = new_figure('Fig6e sweep');
subplot(1,2,1); hold on;
mk = {'o-','s-','d-','^-','v--','x--'};
for v = 1:numel(VAR)
    m = strcmp(rowsR(:,2), VAR{v});
    plot(cell2mat(rowsR(m,1)), cell2mat(rowsR(m,6)), mk{v}, 'LineWidth', 1.2);
end
xline_x = [8 14 20];
for x = xline_x, plot([x x], [0 100], ':', 'Color', [0.6 0.6 0.6]); end
xlabel('boundary separation (px)'); ylabel('both resolved (% rows)');
title('Resolution vs separation'); legend(VAR, 'Location', 'southeast'); grid on;
subplot(1,2,2); hold on;
for v = 1:numel(VAR)
    m = strcmp(rowsW(:,2), VAR{v});
    plot(cell2mat(rowsW(m,1)), cell2mat(rowsW(m,3)), mk{v}, 'LineWidth', 1.2);
end
xlabel('boundary separation (px)'); ylabel('mean detection-run width (px)');
title('Run width before LRT'); grid on;
save_png(f, fullfile(FIG, 'fig6e_sweep_runwidth.png'));
fprintf('\nSaved results_fig6e_extended.csv, run_width_fig6e.csv, fig6e_sweep_runwidth.png\n');
