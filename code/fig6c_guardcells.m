%% ========================================================================
%  fig6c_guardcells.m  —  Task 7: does the CFAR window explain Fig 6c?
%
%  Our Fig 6c reproduction reaches ~86% TPR against the published 100%.
%  os_cfar.m uses the full s x s window, so the cell under test (CUT) and
%  its line neighbours help set their own threshold. The paper places the
%  CUT between two reference kernels. This tests that hypothesis.
%
%  CFAR variants
%    orig  full s x s window (current method)
%    g0    two s x s/2 kernels, CUT column excluded
%    g1    3-column guard strip (covers a 3-px line)
%    g2    5-column guard strip
%  All keep s^2 reference samples, so alpha from calibration.mat is valid.
%  Each variant must pass the same false-alarm gate before it is used.
%
%  Requires ../data/simulated/fig6abcd.mat and calibration.mat
%  (both written by rosed_pipeline.m v2).
% =========================================================================
clear; clc; close all;
import rosed.*          % the shared library, see code\rosed.m

sz = 256;  LRT_N = 14;  guard = 2;  P_FA = 1e-3;
TAB = '../results/tables';
C = load('calibration.mat');
D = load('../data/simulated/fig6abcd.mat');
scales = C.scales;
CF = {'orig', 'g0', 'g1', 'g2'};
G  = [NaN, 0, 1, 2];

%% Gate: false-alarm rate on fresh homogeneous clutter
fprintf('=== False-alarm gate per CFAR variant ===\n');
rng(99);
[~, Pc0] = simulate_speckle(sz, sz, -6);
gate = true(1, numel(CF));
for c = 1:numel(CF)
    for i = 1:numel(scales)
        s = scales(i);  k = round(0.75*s^2);  b = s + 3;
        [dp, dn] = cfar_variant(Pc0, s, k, C.alpha_pos(i), C.alpha_neg(i), G(c));
        rp = nnz(dp(b:end-b, b:end-b)) / numel(dp(b:end-b, b:end-b));
        rn = nnz(dn(b:end-b, b:end-b)) / numel(dn(b:end-b, b:end-b));
        okc = rp >= P_FA/10 && rp <= P_FA*10 && rn >= P_FA/10 && rn <= P_FA*10;
        gate(c) = gate(c) && okc;
        fprintf('%-5s s=%2d  pos %.2e  neg %.2e  %s\n', CF{c}, s, rp, rn, iff(okc, 'ok', 'FAIL'));
    end
end

%% Evaluate
tests = {D.Pb, D.gtB, 'Fig6b'; D.Pc, D.gtC, 'Fig6c'; D.Pa, D.gtA, 'Fig6a'};
VAR = {'fixed-14', 'intersection'};
hdr = {'CFAR','Image','Variant','FN','FP','TPR'};
rows = cell(0, 6);
hdrD = {'CFAR','Variant','BoundaryTPR','BrightTPR','DarkTPR','FP'};
rowsD = cell(0, 6);

for c = 1:numel(CF)
    if ~gate(c)
        fprintf('\n%s failed the gate; skipped.\n', CF{c});
        continue;
    end
    fprintf('\n=== CFAR %s ===\n', CF{c});
    for t = 1:size(tests, 1)
        [S, cf] = stack_with_cfar(tests{t,1}, C, G(c));
        for v = 1:numel(VAR)
            final = fuse_variant(tests{t,1}, S, VAR{v}, LRT_N, guard, cf);
            [FN, FP, TPR] = evaluate_perrow(final, tests{t,2}, 2);
            fprintf('%-6s %-13s FN %4d  FP %5d  TPR %6.1f\n', tests{t,3}, VAR{v}, FN, FP, TPR);
            rows(end+1, :) = {CF{c}, tests{t,3}, VAR{v}, FN, FP, TPR}; %#ok<SAGROW>
        end
    end
    [S, cf] = stack_with_cfar(D.Pd, C, G(c));
    for v = 1:numel(VAR)
        final = fuse_variant(D.Pd, S, VAR{v}, LRT_N, guard, cf);
        Rm = evaluate_multi(final, D.gtD, 2);
        fprintf('Fig6d  %-13s bound %5.1f  bright %5.1f  dark %5.1f  FP %5d\n', VAR{v}, ...
                Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total);
        rowsD(end+1, :) = {CF{c}, VAR{v}, Rm.boundary.TPR, Rm.bright.TPR, Rm.dark.TPR, Rm.FP_total}; %#ok<SAGROW>
    end
end
write_csv(fullfile(TAB, 'cfar_guard_6abc.csv'), hdr, rows);
write_csv(fullfile(TAB, 'cfar_guard_6d.csv'), hdrD, rowsD);

fprintf('\nPublished Fig 6c: FN 0 / FP 641 / TPR 100.\n');
fprintf('If a guard variant lifts Fig6c TPR toward 100 without breaking Fig6b\n');
fprintf('or inflating FP, the full-window CFAR explains the shortfall.\n');
fprintf('Saved cfar_guard_6abc.csv, cfar_guard_6d.csv\n');

%% ---- local helpers ------------------------------------------------------
function [dp, dn] = cfar_variant(P, s, k, ap, an, g)
    import rosed.*          % scoped per function, so repeated here
    if isnan(g)
        dp = os_cfar(P, s, k, ap, 'positive');
        dn = os_cfar(P, s, k, an, 'negative');
    else
        dp = os_cfar_guard(P, s, g, k, ap, 'positive');
        dn = os_cfar_guard(P, s, g, k, an, 'negative');
    end
end

function [S, cf] = stack_with_cfar(P, C, g)
    import rosed.*          % scoped per function, so repeated here
    S = build_scale_stack(P, C.scales, C.taus, C.alpha_pos, C.alpha_neg);
    if ~isnan(g)
        for i = 1:numel(S)
            s = S(i).scale;  k = round(0.75*s^2);
            [S(i).cfar_p, S(i).cfar_n] = cfar_variant(P, s, k, C.alpha_pos(i), C.alpha_neg(i), g);
        end
    end
    cf = cfar_fused(S);
end

function out = iff(c, a, b)
    if c, out = a; else, out = b; end
end
