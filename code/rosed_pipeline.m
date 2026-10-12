%% ========================================================================
%  rosed_pipeline.m  (v2)  —  Simulated reproduction, Fig 6a-d
%  Aditya Rane (B129) / Jay Gandhi (B124)
%  Base: Hashem et al., "Robust SAR Edge Detection", EuRAD 2025
%
%  CHANGES FROM v1
%   - Local functions removed. Every function now lives in its own file in
%     code/, so all experiments call identical code.
%   - Test images and calibration are saved for the other experiments:
%       ../data/simulated/fig6abcd.mat
%       calibration.mat  (here, as before) and ../results/tables/calibration.mat
%   - Results are compared automatically against your saved v1 CSVs.
%     If anything differs, new results go to *_NEW.csv and nothing is
%     overwritten.
%
%  Random numbers are consumed in exactly the same order as v1, so with
%  rng(42) every image and every number must reproduce v1 exactly.
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m
rng(42);

scales = [8 14 20];
P_FA   = 1e-3;
sz     = 256;
LRT_N  = 14;
guard  = 2;
TAB    = '../results/tables';
SIM    = '../data/simulated';
FIG    = '../results/figures';
if ~exist(TAB, 'dir'), mkdir(TAB); end
if ~exist(SIM, 'dir'), mkdir(SIM); end
if ~exist(FIG, 'dir'), mkdir(FIG); end

%% SECTION 1 — Simulator verification
fprintf('=== Section 1: Simulator verification ===\n');
[A0, P0] = simulate_speckle(512, 512, -10);
fprintf('mean power   target %.4f | actual %.4f\n', 10^(-10/10), mean(P0(:)));
fprintf('amplitude CV target %.4f | actual %.4f\n', sqrt(4/pi-1), std(A0(:))/mean(A0(:)));
fprintf('intensity CV target %.4f | actual %.4f\n', 1.0, std(P0(:))/mean(P0(:)));

%% SECTION 2 — Test images
fprintf('\n=== Section 2: Test images ===\n');
[Pa, gtA] = make_fig6a(sz);
[Pb, gtB] = make_fig6b(sz);
[Pc, gtC] = make_fig6c(sz);

f = new_figure('Test images');
subplot(1,3,1); imagesc(sqrt(Pa)); axis image; title('Fig 6a — clutter boundary');
subplot(1,3,2); imagesc(sqrt(Pb)); axis image; title('Fig 6b — bright target');
subplot(1,3,3); imagesc(sqrt(Pc)); axis image; title('Fig 6c — dark target');
colormap gray;
save_png(f, fullfile(FIG, 'test_images_6abc.png'));

%% SECTION 3 — ROM threshold calibration
fprintf('\n=== Section 3: ROM calibration ===\n');
taus = zeros(1, numel(scales));
for i = 1:numel(scales)
    taus(i) = calibrate_rom(scales(i), P_FA, 3, sz);
end

%% SECTION 3b — OS-CFAR alpha
alpha_pos = zeros(1, numel(scales));
alpha_neg = zeros(1, numel(scales));
for i = 1:numel(scales)
    N = scales(i)^2;
    k = round(0.75*N);
    alpha_pos(i) = solve_alpha(N, k, P_FA);
    alpha_neg(i) = solve_alpha(N, k, 1 - P_FA);
end
fprintf('\n%6s %10s %10s %12s\n', 'scale', 'tau', 'alpha_pos', 'alpha_neg');
for i = 1:numel(scales)
    fprintf('%6d %10.4f %10.4f %12.4e\n', scales(i), taus(i), alpha_pos(i), alpha_neg(i));
end
save('calibration.mat', 'taus', 'scales', 'P_FA', 'alpha_pos', 'alpha_neg', '-v7');
save(fullfile(TAB, 'calibration.mat'), 'taus', 'scales', 'P_FA', 'alpha_pos', 'alpha_neg', '-v7');

%% SECTION 4 — GATE: CFAR false-alarm rate
fprintf('\n=== Section 4: CFAR calibration check ===\n');
gate_pass = true;
for i = 1:numel(scales)
    [rp, rn] = verify_cfar(scales(i), alpha_pos(i), alpha_neg(i), sz);
    fprintf('scale %2d: pos=%.2e  neg=%.2e  (target %.1e)\n', scales(i), rp, rn, P_FA);
    if rp < P_FA/10 || rp > P_FA*10
        gate_pass = false;
    end
end
if ~gate_pass
    error('CFAR calibration FAILED. Do not interpret results below.');
end
fprintf('GATE PASSED.\n');

%% SECTION 5 — Reproduction, Fig 6a/b/c
fprintf('\n=== Section 5: Reproduction (Fig 6a/b/c) ===\n');
tests  = {Pa, gtA, 'Fig6a'; Pb, gtB, 'Fig6b'; Pc, gtC, 'Fig6c'};
vnames = {'fixed-14', 'union', 'intersection', 'corroborated'};
hdr5   = {'Image', 'Variant', 'FN', 'FP', 'TPR'};
rows5  = cell(0, 5);

for t = 1:3
    Pimg = tests{t,1};  gtcols = tests{t,2};  name = tests{t,3};
    S    = build_scale_stack(Pimg, scales, taus, alpha_pos, alpha_neg);
    cf   = cfar_fused(S);
    fprintf('\n--- %s ---\n%-14s %8s %10s %8s\n', name, 'variant', 'FN', 'FP', 'TPR%');
    for v = 1:4
        [final, ~, sel] = fuse_variant(Pimg, S, vnames{v}, LRT_N, guard, cf);
        [FN, FP, TPR] = evaluate_perrow(final, gtcols, 2);
        fprintf('%-14s %8d %10d %8.1f\n', vnames{v}, FN, FP, TPR);
        rows5(end+1, :) = {name, vnames{v}, FN, FP, TPR}; %#ok<SAGROW>
        if v == 4
            f = new_figure(['Scale selection — ' name]);
            imagesc(sel); colorbar; set(gca, 'CLim', [0 3]);
            title(['Scale selection — ' name ' (0=none, 1=fine, 2=mid, 3=coarse)']);
            save_png(f, fullfile(FIG, ['scale_selection_' name '.png']));
        end
    end
end

%% SECTION 5b — Mixed scene, Fig 6d
fprintf('\n=== Section 5b: Mixed scene Fig6d ===\n');
[Pd, gtD] = make_fig6d(sz);
Sd = build_scale_stack(Pd, scales, taus, alpha_pos, alpha_neg);
cfd = cfar_fused(Sd);
hdr5b = {'Image', 'Variant', 'BoundaryTPR', 'BrightTPR', 'DarkTPR', 'FP'};
rows5b = cell(0, 6);
fprintf('%-14s %11s %11s %11s %10s\n', 'variant', 'bound TPR', 'bright TPR', 'dark TPR', 'FP');
for v = 1:4
    [final, ~, sel] = fuse_variant(Pd, Sd, vnames{v}, LRT_N, guard, cfd);
    Rm = evaluate_multi(final, gtD, 2);
    fprintf('%-14s %11.1f %11.1f %11.1f %10d\n', vnames{v}, ...
            Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total);
    rows5b(end+1, :) = {'Fig6d', vnames{v}, Rm.boundary.TPR, Rm.bright.TPR, ...
                        Rm.dark.TPR, Rm.FP_total}; %#ok<SAGROW>
    if v == 4
        f = new_figure('Scale selection — Fig6d');
        imagesc(sel); colorbar; set(gca, 'CLim', [0 3]);
        title('Scale selection — Fig6d (0=none, 1=fine, 2=mid, 3=coarse)');
        save_png(f, fullfile(FIG, 'scale_selection_Fig6d.png'));
    end
end

save(fullfile(SIM, 'fig6abcd.mat'), 'Pa', 'gtA', 'Pb', 'gtB', 'Pc', 'gtC', ...
     'Pd', 'gtD', 'sz', '-v7');
fprintf('\nSaved test images to %s/fig6abcd.mat\n', SIM);

%% SECTION 6 — Classical baseline, unthinned Sobel (v1 method, kept for regression)
fprintf('\n=== Section 6: Unthinned Sobel (v1) ===\n');
allT = [tests; {Pd, gtD, 'Fig6d'}];
for t = 1:4
    e = sobel_edge(norm01(sqrt(allT{t,1})));
    if isstruct(allT{t,2})
        Rm = evaluate_multi(e, allT{t,2}, 2);
        fprintf('%-6s bound=%.1f%% bright=%.1f%% dark=%.1f%% FP=%d\n', allT{t,3}, ...
                Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total);
    else
        [FN, FP, TPR] = evaluate_perrow(e, allT{t,2}, 2);
        fprintf('%-6s FN=%d FP=%d TPR=%.1f%%\n', allT{t,3}, FN, FP, TPR);
    end
end
fprintf('Thinned Sobel, LoG and Canny: run classical_baselines.m\n');

%% SECTION 7 — Regression against v1 results, then save
fprintf('\n=== Section 7: Regression check against v1 ===\n');
ok1 = compare_to_baseline({fullfile(TAB,'results_table.csv'), 'results_table.csv'}, ...
                          hdr5, rows5, {'Image','Variant'}, {'FN','FP','TPR'}, 1e-3);
ok2 = compare_to_baseline({fullfile(TAB,'results_fig6d.csv'), 'results_fig6d.csv'}, ...
                          hdr5b, rows5b, {'Image','Variant'}, ...
                          {'BoundaryTPR','BrightTPR','DarkTPR','FP'}, 1e-3);
REGRESSION_OK = ok1 && ok2;
if REGRESSION_OK
    write_csv(fullfile(TAB, 'results_table.csv'), hdr5, rows5);
    write_csv(fullfile(TAB, 'results_fig6d.csv'), hdr5b, rows5b);
    fprintf('\nREGRESSION PASSED. Results saved to %s\n', TAB);
else
    write_csv(fullfile(TAB, 'results_table_NEW.csv'), hdr5, rows5);
    write_csv(fullfile(TAB, 'results_fig6d_NEW.csv'), hdr5b, rows5b);
    fprintf('\nREGRESSION FAILED. v1 files untouched; new results in *_NEW.csv.\n');
    fprintf('Do not continue until this is understood.\n');
end
