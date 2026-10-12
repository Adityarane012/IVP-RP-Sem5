function diagnose_patches(cropdir, pattern, s, P_FA)
% DIAGNOSE_PATCHES  Per-patch tau, to test whether pooling is valid.
%
% Pooling r-hat across patches assumes they share a distribution. If the
% per-patch taus scatter widely, that assumption fails and the pooled
% percentile is set by whichever few patches have the heaviest tails
% rather than by the population.
%
% Also reports the fraction of very small pixel values in each patch.
% These uint8 products bottom out near 1, so a half-window median can
% land on a small integer; a median of 2 against a median of 25 yields
% r-hat = 0.08 from quantization alone, with no edge present.
%
% Usage:
%   diagnose_patches('../data/crops', '*clean*.tif', 14, 1e-3);
    import rosed.*          % the shared library, see code\rosed.m

    files = dir(fullfile(cropdir, pattern));
    if isempty(files)
        error('No files matching %s in %s', pattern, cropdir);
    end

    n = numel(files);
    tau  = zeros(1,n);
    cv   = zeros(1,n);
    mu   = zeros(1,n);
    lowf = zeros(1,n);
    p01  = zeros(1,n);

    fprintf('\n=== PER-PATCH DIAGNOSTIC  (scale %d, P_FA %.0e) ===\n\n', s, P_FA);
    fprintf('%-42s %6s %6s %7s %8s %8s\n', ...
            'patch', 'mean', 'CV', 'tau', 'p<=3%%', 'rhat_p1');

    b = s;
    for i = 1:n
        f = fullfile(files(i).folder, files(i).name);
        A = double(imread(f));

        mu(i)   = mean(A(:));
        cv(i)   = std(A(:))/mu(i);
        lowf(i) = 100 * nnz(A <= 3) / numel(A);

        rh = rom_map(A, s);
        v  = rh(b:end-b, b:end-b);
        v  = sort(v(:));
        tau(i) = v(max(1, round(P_FA*numel(v))));
        p01(i) = v(max(1, round(0.01*numel(v))));   % 1st percentile

        [~, nm] = fileparts(files(i).name);
        fprintf('%-42s %6.1f %6.3f %7.4f %8.2f %8.4f\n', ...
                nm(max(1,end-41):end), mu(i), cv(i), tau(i), lowf(i), p01(i));
    end

    fprintf('\n--- SPREAD ---\n');
    fprintf('tau   : min %.4f  median %.4f  max %.4f   ratio max/min = %.1f\n', ...
            min(tau), median(tau), max(tau), max(tau)/max(min(tau),eps));
    fprintf('CV    : min %.3f  median %.3f  max %.3f\n', ...
            min(cv), median(cv), max(cv));
    fprintf('mean  : min %.1f  median %.1f  max %.1f\n', ...
            min(mu), median(mu), max(mu));
    fprintf('%% pix <= 3 : min %.2f  median %.2f  max %.2f\n', ...
            min(lowf), median(lowf), max(lowf));

    fprintf('\n--- CORRELATIONS ---\n');
    fprintf('tau vs patch mean : r = %+.3f\n', corr2(mu, tau));
    fprintf('tau vs patch CV   : r = %+.3f\n', corr2(cv, tau));
    fprintf('tau vs %% low pix  : r = %+.3f\n', corr2(lowf, tau));

    fprintf('\n--- READING THIS ---\n');
    fprintf('tau max/min under ~3    : patches comparable, pooling is valid\n');
    fprintf('tau max/min over ~10    : patches differ, pooled tail is driven\n');
    fprintf('                          by a few patches, not the population\n');
    fprintf('strong -ve corr with mean, or +ve with %% low pix:\n');
    fprintf('                          quantization is setting the tail\n');

    figure('Name','Per-patch tau','Color','w','Position',[100 100 1100 380]);
    subplot(1,3,1); scatter(mu, tau, 40, 'filled');
        xlabel('patch mean'); ylabel('tau'); grid on; title('tau vs mean');
    subplot(1,3,2); scatter(cv, tau, 40, 'filled');
        xlabel('patch CV'); ylabel('tau'); grid on; title('tau vs CV');
    subplot(1,3,3); scatter(lowf, tau, 40, 'filled');
        xlabel('% pixels <= 3'); ylabel('tau'); grid on; title('tau vs low-value fraction');
end

function r = corr2(x, y)
    x = x(:) - mean(x); y = y(:) - mean(y);
    r = (x'*y) / max(sqrt((x'*x)*(y'*y)), eps);
end
