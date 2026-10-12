%% ========================================================================
%  real_calibrate_neg.m  —  Task 2: negative CFAR constant on real data
%
%  The real-data calibration (run_pooled_calibration.m) found only the
%  positive CFAR constant. The negative detector needs its own.
%
%  TIES: uint8 data make A/ghat take repeated exact values, so the rate
%  cannot always equal P_FA. Among achievable thresholds, the one whose
%  actual rate is closest to P_FA (in log terms) is chosen and reported.
%
%  EXACT METHOD, NO SEARCH
%  The negative test flags a pixel when  A < alpha * ghat,  which is the
%  same as  A/ghat < alpha.  So on homogeneous clutter the false-alarm rate
%  equals P_FA exactly when alpha is the P_FA-quantile of the pooled ratio
%  A/ghat. Likewise alpha_pos is the (1 - P_FA)-quantile. This gives both
%  constants directly, and reproducing the bisected alpha_pos (about 1.97)
%  is a built-in check that the method is consistent.
%
%  Operates on AMPLITUDE, like the rest of the real-data pipeline.
%  Uses the 12-region pool in ../data/crops/*clean*.tif.
%  Appends alphas_neg_real to calibration_real_pooled.mat.
% =========================================================================
clear; clc;
import rosed.*          % the shared library, see code\rosed.m

CALF = '../results/tables/calibration_real_pooled.mat';
R = load(CALF);
SCALES = R.SCALES;  P_FA = R.P_FA;
files = dir(fullfile('../data/crops', '*clean*.tif'));
if isempty(files), error('No *clean*.tif regions in ../data/crops'); end
fprintf('regions: %d\n', numel(files));

alphas_neg_real = zeros(1, numel(SCALES));
alpha_pos_check = zeros(1, numel(SCALES));
rate_neg = zeros(1, numel(SCALES));

for si = 1:numel(SCALES)
    s = SCALES(si);  k = round(0.75*s^2);  b = s;
    ratio = [];
    for i = 1:numel(files)
        A = double(imread(fullfile(files(i).folder, files(i).name)));
        g = ordstat_filter(A, k, true(s, s));        % same reference as os_cfar.m
        q = A(b:end-b, b:end-b) ./ max(g(b:end-b, b:end-b), eps);
        ratio = [ratio; q(:)]; %#ok<AGROW>
    end
    [aN, rN] = tie_aware(ratio, P_FA, 'below');
    [aP, rP] = tie_aware(ratio, P_FA, 'above');
    alphas_neg_real(si) = aN;  rate_neg(si) = rN;  alpha_pos_check(si) = aP;
    fprintf('s=%2d  alpha_neg %.5f  achieved rate %.2e | alpha_pos check %.3f (rate %.2e) vs saved %.3f\n', ...
            s, aN, rN, aP, rP, R.alphas_real(si));
    if rN < P_FA/3 || rN > 3*P_FA
        fprintf('   WARNING: no threshold gives a rate near %.0e. 8-bit quantisation makes\n', P_FA);
        fprintf('   A/ghat discrete near zero. Record the achieved rate in the paper.\n');
    end
end

% Rayleigh-amplitude reference
a_theory = sqrt(-log(1 - P_FA)) / sqrt(log(4));
fprintf('\nRayleigh amplitude prediction: alpha_neg = %.4f (asymptotic)\n', a_theory);

R.alphas_neg_real = alphas_neg_real;
R.rate_neg        = rate_neg;
R.alpha_pos_check = alpha_pos_check;
save(CALF, '-struct', 'R', '-v7');
fprintf('Saved alphas_neg_real to %s\n', CALF);

function [alpha, rate] = tie_aware(ratio, P_FA, side)
% Choose the threshold whose ACTUAL strict rate is closest to P_FA.
%   below:  rate = mean(ratio <  alpha)   (negative detector)
%   above:  rate = mean(ratio >  alpha)   (positive detector)
    q = sort(ratio(:));  n = numel(q);
    if strcmp(side, 'below')
        [u, ia] = unique(q, 'first');
        rates = (ia - 1) / n;              % count strictly below each value
    else
        [u, ib] = unique(q, 'last');
        rates = (n - ib) / n;              % count strictly above each value
    end
    ok = rates > 0;
    if ~any(ok)
        alpha = u(end);  rate = 0;  return;
    end
    err = abs(log(rates(ok) / P_FA));
    uu = u(ok);  rr = rates(ok);
    [~, j] = min(err);
    alpha = uu(j);  rate = rr(j);
end
