classdef rosed
% ROSED  The shared library for this project, as one file.
%
%   Multi-Scale Intersection for Improved Edge Separation in
%   Robust SAR Edge Detection.  Jay Gandhi, Aditya Rane.
%
%   Every function that more than one experiment uses lives here, so
%   that all results come from one copy of the code. The experiments
%   themselves stay as separate scripts, because each one is a
%   separate piece of work.
%
%   To use it, put one line at the top of the caller:
%
%       import rosed.*
%
%   after which every function below is called by its plain name,
%   exactly as when each one was its own file. A qualified call such
%   as rosed.rom_map(P, 14) also works, including inside parfor.
%
%   Chapters:
%     1  Speckle simulation and the simulated test images         6 functions
%     2  The detector                                             9 functions
%     3  Combining the three window sizes                         3 functions
%     4  Threshold calibration and verification                   4 functions
%     5  Evaluation metrics                                       4 functions
%     6  Regression gate and file input and output                3 functions
%     7  Small utilities                                          6 functions
%
%   Assembled from the 35 separate files by
%   the project's build script. The originals are kept in
%   the project's local working copy.

    methods (Static)

    % ====================================================================
    %  CHAPTER 1  Speckle simulation and the simulated test images
    %  Follows the physical model in [1]: I and Q are Gaussian, so the
    %  intensity is exponential and the amplitude Rayleigh.
    % ====================================================================

        function [A, P] = simulate_speckle(rows, cols, dB)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            sigma = 10^(dB/20);
            I = randn(rows, cols) * sigma/sqrt(2);
            Q = randn(rows, cols) * sigma/sqrt(2);
            P = I.^2 + Q.^2;    % exponential, mean sigma^2
            A = sqrt(P);        % Rayleigh, display only
        end

        function [P, gtcols] = make_fig6a(sz)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % The paper states mean powers of -10 dB and -8 dB AND a contrast ratio
        % rho = 10 in the amplitude image. These are inconsistent: a 2 dB power
        % difference gives a power ratio of 1.585 and therefore an amplitude
        % ratio of 1.26, not 10. Obtaining rho = 10 requires a 20 dB difference.
        % We adopt rho = 10, since only that reproduces the published 100% TPR.
        % Documented as a reproducibility finding.
            import rosed.*
            [~, PL] = simulate_speckle(sz, sz/2, -10);
            [~, PR] = simulate_speckle(sz, sz/2, +10);
            P = [PL, PR];
            gtcols = [sz/2, sz/2+1];
        end

        function [P, gtcols] = make_fig6b(sz)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            import rosed.*
            [~, P] = simulate_speckle(sz, sz, -6);
            c = round(sz/2);
            [~, tgt] = simulate_speckle(sz, 3, 4);
            P(:, c-1:c+1) = tgt;
            gtcols = c-1:c+1;
        end

        function [P, gtcols] = make_fig6c(sz)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            import rosed.*
            [~, P] = simulate_speckle(sz, sz, -6);
            c = round(sz/2);
            [~, tgt] = simulate_speckle(sz, 3, -16);
            P(:, c-1:c+1) = tgt;
            gtcols = c-1:c+1;
        end

        function [P, gtcols] = make_fig6d(sz)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Broad clutter boundary AND narrow line targets in a single scene.
        % No such image exists in the paper's evaluation set. A fixed kernel
        % must compromise here, so this is the only regime in which adaptive
        % scale selection can demonstrate a benefit.
            import rosed.*
            [~, PL] = simulate_speckle(sz, sz/2, -10);
            [~, PR] = simulate_speckle(sz, sz/2, +10);
            P = [PL, PR];

            [~, t1] = simulate_speckle(sz, 3, 4);     % bright line, dark half
            P(:, 60:62) = t1;

            [~, t2] = simulate_speckle(sz, 3, -6);    % dark line, bright half
            P(:, 190:192) = t2;

            gtcols = struct('boundary', [sz/2, sz/2+1], ...
                            'bright',   60:62, ...
                            'dark',     190:192);
        end

        function [P, gt] = make_fig6e(sz, gap)
        % [shared library] Extracted unchanged from fig6e_experiment.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Two clutter boundaries separated by `gap` pixels.
        % Dark background, a brighter band of width `gap`, dark again. Contrast
        % uses the same rho = 10 amplitude ratio adopted for Fig6a, i.e. a 20 dB
        % power difference, since the paper's stated -10/-8 dB values are
        % inconsistent with its stated rho.
            import rosed.*
            c  = round(sz/2);
            x1 = c - floor(gap/2);          % first boundary
            x2 = x1 + gap;                  % second boundary

            [~, PL] = simulate_speckle(sz, x1,            -10);
            [~, PM] = simulate_speckle(sz, gap,           +10);
            [~, PR] = simulate_speckle(sz, sz - x1 - gap, -10);
            P = [PL, PM, PR];

            gt = struct('edge1', [x1, x1+1], 'edge2', [x2, x2+1]);
        end

    % ====================================================================
    %  CHAPTER 2  The detector
    %  ROM for clutter boundaries, OS-CFAR for bright and dark lines, and
    %  the LRT to place one edge per detection run.
    % ====================================================================

        function rhat = rom_map(P, s)
        % Ratio of Medians. Window split into left and right halves along x.
        % r = median(left)/median(right);  rhat = min(r, 1/r) in [0,1].
        % Low rhat means strong edge evidence.
            import rosed.*
            h = s/2;
            domL = false(s,s);  domL(:, 1:h)     = true;
            domR = false(s,s);  domR(:, h+1:end) = true;
            ord  = round(s*h/2);

            m1 = ordstat_filter(P, ord, domL);
            m2 = ordstat_filter(P, ord, domR);

            r = m1 ./ max(m2, eps);
            rhat = min(r, 1./max(r, eps));
            rhat(~isfinite(rhat)) = 1;
        end

        function det = os_cfar(P, s, k, alpha, mode)
        % Order-Statistics CFAR. ghat is the k-th order statistic of the window;
        % a pixel is flagged when it exceeds (positive) or falls below (negative)
        % alpha times that reference.
            import rosed.*
            dom  = true(s,s);
            ghat = ordstat_filter(P, k, dom);
            if strcmp(mode, 'positive')
                det = P > alpha .* ghat;
            else
                det = P < alpha .* ghat;
            end
        end

        function det = os_cfar_guard(P, s, g, k, alpha, mode)
        % OS_CFAR_GUARD  OS-CFAR with the cell under test outside its own reference.
        %
        % The base paper places the cell under test (CUT) between two reference
        % kernels. os_cfar.m instead uses the full s x s window, so the CUT and its
        % immediate neighbours help set their own threshold. For a 3-px dark line
        % that lowers the reference and suppresses detection.
        %
        % Reference domain here: two s x (s/2) kernels separated by a guard strip of
        % 2g+1 columns centred on the CUT. The reference still holds s^2 samples, so
        % the CFAR constant alpha solved for (N = s^2, k) remains valid.
        %
        %   g = 0  excludes only the CUT column
        %   g = 1  excludes a 3-column strip (covers a 3-px line)
        %   g = 2  excludes a 5-column strip
            import rosed.*
            h   = s/2;
            dom = [true(s, h), false(s, 2*g + 1), true(s, h)];
            ghat = ordstat_filter(P, k, dom);
            if strcmp(mode, 'positive')
                det = P > alpha .* ghat;
            else
                det = P < alpha .* ghat;
            end
        end

        function out = ordstat_filter(P, order, dom)
        % Uses the compiled ordfilt2. Verified bit-identical to the hand-written
        % version by test_ordfilt2.m across all nine domains the pipeline uses.
            out = ordfilt2(P, order, dom, 'symmetric');
        end

        function best = lrt_1d(profile, guard)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Log-likelihood ratio for a step in exponential data.
            N    = numel(profile);
            Stot = sum(profile);
            bestT = -Inf;
            best  = round(N/2);
            for i = (guard+1):(N-guard-1)
                S1 = sum(profile(1:i));
                S2 = Stot - S1;
                n1 = i;
                n2 = N - i;
                if S1 <= 0 || S2 <= 0, continue; end
                T = N*log(Stot/N) - n1*log(S1/n1) - n2*log(S2/n2);
                if T > bestT
                    bestT = T;
                    best  = i;
                end
            end
        end

        function loc = localize_rom(P, rom_edge, N, guard)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % LRT localizes clutter boundaries only. ROM gives thick detections;
        % each thick run collapses to its single most likely position.
            import rosed.*
            [R, C] = size(P);
            loc = false(R, C);
            for y = 1:R
                cols = find(rom_edge(y,:));
                if isempty(cols), continue; end
                cols = cols(:)';
                brk = [0, find(diff(cols) > 2), numel(cols)];
                for c = 1:numel(brk)-1
                    seg = cols(brk(c)+1 : brk(c+1));
                    x0  = round(mean(seg));
                    lo  = max(1, x0 - floor(N/2));
                    hi  = min(C, lo + N - 1);
                    lo  = max(1, hi - N + 1);
                    if hi - lo + 1 < 2*guard + 3, continue; end
                    i_star = lrt_1d(P(y, lo:hi), guard);
                    loc(y, lo + i_star - 1) = true;
                end
            end
        end

        function S = build_scale_stack(P, scales, taus, alpha_pos, alpha_neg)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            import rosed.*
            S = struct('scale', {}, 'rhat', {}, 'rom', {}, 'cfar_p', {}, 'cfar_n', {});
            for i = 1:numel(scales)
                s = scales(i);
                N = s*s;
                k = round(0.75*N);
                S(i).scale  = s;
                S(i).rhat   = rom_map(P, s);
                S(i).rom    = S(i).rhat < taus(i);
                S(i).cfar_p = os_cfar(P, s, k, alpha_pos(i), 'positive');
                S(i).cfar_n = os_cfar(P, s, k, alpha_neg(i), 'negative');
            end
        end

        function C = cfar_fused(S)
        % CFAR_FUSED  OR of positive and negative CFAR maps across all scales.
            C = false(size(S(1).rom));
            for i = 1:numel(S)
                C = C | S(i).cfar_p | S(i).cfar_n;
            end
        end

        function [final, romsel, sel] = fuse_variant(P, S, name, LRT_N, guard, cfar_all)
        % FUSE_VARIANT  Final ROSED edge map for one scale-combination variant.
        %
        %   P        image used by the LRT (intensity)
        %   S        scale stack from build_scale_stack
        %   name     'fixed-N' (e.g. 'fixed-14') | 'union' | 'intersection' | 'corroborated'
        %   cfar_all CFAR map fused across scales. Identical for every variant, so
        %            differences between variants come from ROM selection alone.
        %
        % Matches the per-variant logic used in rosed_pipeline.m and
        % fig6e_experiment.m exactly.
            import rosed.*
            switch name
                case 'union'
                    [romsel, sel] = select_union(S);
                case 'intersection'
                    [romsel, sel] = select_intersection(S);
                case 'corroborated'
                    [romsel, sel] = select_corroborated(S, 2);
                otherwise
                    % 'fixed-N' : single scale N from the stack
                    tok = regexp(name, '^fixed-(\d+)$', 'tokens', 'once');
                    if isempty(tok), error('Unknown variant %s', name); end
                    sN = str2double(tok{1});
                    iN = find([S.scale] == sN, 1);
                    if isempty(iN), error('%s requires scale %d in the stack', name, sN); end
                    romsel = S(iN).rom;
                    sel    = iN * double(romsel);
                    if sN == 14, sel = 2 * double(romsel); end   % v1 convention
            end
            final = localize_rom(P, romsel, LRT_N, guard) | cfar_all;
        end

    % ====================================================================
    %  CHAPTER 3  Combining the three window sizes
    %  select_intersection is the contribution of this work: the logical AND
    %  of the ROM maps, applied before the run grouping and the LRT.
    % ====================================================================

        function [E, sel] = select_intersection(S)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            E = S(1).rom;
            for i = 2:numel(S)
                E = E & S(i).rom;
            end
            sel = numel(S) * double(E);
        end

        function [E, sel] = select_union(S)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            E   = false(size(S(1).rom));
            sel = zeros(size(E));
            for i = numel(S):-1:1
                newly = S(i).rom & ~E;
                sel(newly) = i;
                E = E | S(i).rom;
            end
        end

        function [E, sel] = select_corroborated(S, radius)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Coarsest scale trusted directly: most samples, most stable median.
        % Finer scales accepted only where CFAR independently agrees, since a
        % genuine narrow structure is both a distributional discontinuity and an
        % intensity outlier, while speckle fooling an undersized median is
        % unlikely to trip both.
            import rosed.*
            nS  = numel(S);
            E   = S(nS).rom;
            sel = zeros(size(E));
            sel(E) = nS;
            for i = nS-1:-1:1
                corr   = dilate_disk(S(i).cfar_p | S(i).cfar_n, radius);
                accept = S(i).rom & corr & ~E;
                E(accept)   = true;
                sel(accept) = i;
            end
        end

    % ====================================================================
    %  CHAPTER 4  Threshold calibration and verification
    %  Monte Carlo thresholds for the simulated data and pooled thresholds
    %  for the real crops.
    % ====================================================================

        function tau = calibrate_rom(s, P_FA, n_trials, sz)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Monte Carlo null calibration on homogeneous clutter. The null
        % distribution of rhat is tighter at large s, so one shared threshold
        % would silently run the small kernel at a higher false-alarm rate.
            import rosed.*
            r_null = [];
            for t = 1:n_trials
                [~, P] = simulate_speckle(sz, sz, -6);
                r = rom_map(P, s);
                b = s;
                v = r(b:end-b, b:end-b);
                r_null = [r_null; v(:)]; %#ok<AGROW>
            end
            srt = sort(r_null);
            tau = srt(max(1, round(P_FA*numel(srt))));
            fprintf('scale %2d: tau_ROM = %.4f (n=%d)\n', s, tau, numel(r_null));
        end

        function alpha = solve_alpha(N, k, P_FA)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % P_FA = prod_{i=0}^{k-1} (N-i)/(N-i+alpha), exponential clutter.
        % k is a PERCENTILE. The published k = 147 at N = 196 is the 75th
        % percentile, so k must scale with N or the CFAR property is lost.
            f = @(a) prod( (N-(0:k-1)) ./ (N-(0:k-1)+a) ) - P_FA;
            alpha = fzero(f, [1e-8, 1e8]);
        end

        function [rate_p, rate_n] = verify_cfar(s, a_pos, a_neg, sz)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Pure clutter, no targets. Detection rate must approximate P_FA.
        % A rate of exactly zero means the wrong domain is being fed in.
            import rosed.*
            [~, P] = simulate_speckle(sz, sz, -6);
            N = s*s;
            k = round(0.75*N);
            dp = os_cfar(P, s, k, a_pos, 'positive');
            dn = os_cfar(P, s, k, a_neg, 'negative');
            b = s;
            dp = dp(b:end-b, b:end-b);
            dn = dn(b:end-b, b:end-b);
            rate_p = nnz(dp)/numel(dp);
            rate_n = nnz(dn)/numel(dn);
        end

        function [taus, alphas] = pooled_calibrate(cropdir, pattern, scales, P_FA)
        % POOLED_CALIBRATE  Threshold calibration pooled over many clean patches.
        %
        % A single small patch is homogeneous but sample-starved: a 128x128 crop
        % yields only ~14 samples below the 1e-3 percentile, so the estimate is
        % driven by a handful of extreme pixels. A single large crop has enough
        % samples but contains real edges, which contaminate the low tail and
        % push the threshold far too low.
        %
        % Pooling many small clean patches gives both: the sample count of a
        % large region with the homogeneity of a small one.
        %
        % Operates on AMPLITUDE (uint8 values as-is), not intensity.
        %
        % Usage:
        %   [t,a] = pooled_calibrate('../data/crops', 'clean*.tif', [8 14 20], 1e-3);
            import rosed.*

            files = dir(fullfile(cropdir, pattern));
            if isempty(files)
                error('No files matching %s in %s', pattern, cropdir);
            end

            fprintf('\n=== POOLED REAL-DATA CALIBRATION ===\n');
            fprintf('patches : %d\n\n', numel(files));

            % report per-patch statistics first
            allCV = zeros(1, numel(files));
            for i = 1:numel(files)
                A = double(imread(fullfile(files(i).folder, files(i).name)));
                allCV(i) = std(A(:))/mean(A(:));
            end
            fprintf('patch CV: min %.3f  median %.3f  max %.3f\n', ...
                    min(allCV), median(allCV), max(allCV));
            fprintf('(single-look Rayleigh amplitude = 0.523)\n\n');

            taus   = zeros(1, numel(scales));
            alphas = zeros(1, numel(scales));

            for si = 1:numel(scales)
                s = scales(si);
                b = s;

                % ---- pool ROM null samples across every patch ----------------
                pool = [];
                for i = 1:numel(files)
                    A  = double(imread(fullfile(files(i).folder, files(i).name)));
                    rh = rom_map(A, s);
                    v  = rh(b:end-b, b:end-b);
                    pool = [pool; v(:)]; %#ok<AGROW>
                end
                pool = sort(pool);
                idx  = max(1, round(P_FA*numel(pool)));
                taus(si) = pool(idx);

                % ---- pool CFAR detections to find alpha ----------------------
                N = s*s;
                k = round(0.75*N);
                lo = 1; hi = 500;
                for it = 1:30
                    a = (lo + hi)/2;
                    hits = 0; tot = 0;
                    for i = 1:numel(files)
                        A = double(imread(fullfile(files(i).folder, files(i).name)));
                        d = os_cfar(A, s, k, a, 'positive');
                        d = d(b:end-b, b:end-b);
                        hits = hits + nnz(d);
                        tot  = tot  + numel(d);
                    end
                    if hits/tot > P_FA, lo = a; else, hi = a; end
                end
                alphas(si) = (lo + hi)/2;

                fprintf('scale %2d : tau = %.4f  alpha = %.3f  (pooled n = %d, tail idx = %d)\n', ...
                        s, taus(si), alphas(si), numel(pool), idx);
            end

            fprintf('\nSimulated reference: tau 0.2887 / 0.4956 / 0.6263\n');
            fprintf('                     alpha 5.510 / 5.150 / 5.064\n');
            fprintf('\nSANITY CHECK: tau must increase with scale. If it does not,\n');
            fprintf('the patches are not homogeneous enough.\n');
        end

    % ====================================================================
    %  CHAPTER 5  Evaluation metrics
    %  Scoring is per image row, which is what avoids the 50 per cent
    %  ceiling described in the paper.
    % ====================================================================

        function [FN, FP, TPR] = evaluate_perrow(det, gtcols, tol)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Detection rate per row, not pixel overlap. A row counts as detected if
        % any detection falls within tol columns of the true edge. This removes
        % the artificial 50% ceiling from one-pixel-per-row localization.
            [R, C] = size(det);
            lo = max(1, min(gtcols) - tol);
            hi = min(C, max(gtcols) + tol);

            hit     = any(det(:, lo:hi), 2);
            TP_rows = nnz(hit);
            FN      = R - TP_rows;
            TPR     = 100 * TP_rows / R;

            outside = det;
            outside(:, lo:hi) = false;
            FP = nnz(outside);
        end

        function R = evaluate_multi(det, gtstruct, tol)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Fig6d contains three structures at different scales. Scoring them
        % together would average away exactly the effect under test, so each is
        % scored separately. FP counts detections outside every structure.
            [nR, nC] = size(det);
            fields = fieldnames(gtstruct);
            R = struct();
            covered = false(1, nC);

            for f = 1:numel(fields)
                cols = gtstruct.(fields{f});
                lo = max(1, min(cols) - tol);
                hi = min(nC, max(cols) + tol);
                covered(lo:hi) = true;

                hit = any(det(:, lo:hi), 2);
                R.(fields{f}).TPR = 100 * nnz(hit) / nR;
                R.(fields{f}).FN  = nR - nnz(hit);
            end

            outside = det;
            outside(:, covered) = false;
            R.FP_total = nnz(outside);
        end

        function R = evaluate_pair(det, gt, tol)
        % [shared library] Extracted unchanged from fig6e_experiment.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Score the two boundaries separately, and record how often BOTH were
        % resolved in the same row. A detector that merges the pair into one
        % thick response scores well on the individual edges but poorly on
        % "both", which is the quantity that distinguishes the variants.
            [nR, nC] = size(det);

            lo1 = max(1, min(gt.edge1) - tol);  hi1 = min(nC, max(gt.edge1) + tol);
            lo2 = max(1, min(gt.edge2) - tol);  hi2 = min(nC, max(gt.edge2) + tol);

            h1 = any(det(:, lo1:hi1), 2);
            h2 = any(det(:, lo2:hi2), 2);

            R.e1   = 100 * nnz(h1) / nR;
            R.e2   = 100 * nnz(h2) / nR;
            R.both = 100 * nnz(h1 & h2) / nR;

            outside = det;
            outside(:, lo1:hi1) = false;
            outside(:, lo2:hi2) = false;
            R.FP = nnz(outside);
        end

        function R = run_stats(rom_edge, lo, hi)
        % RUN_STATS  Detection-run statistics inside columns lo:hi.
        %
        % This measures the quantity the run-width mechanism is about. localize_rom
        % groups each row's detections into runs (a gap of more than 2 px starts a
        % new run) and emits one edge per run. Two boundaries can only be reported
        % separately if they fall in two separate runs.
        %
        % Returns, over rows that contain at least one run inside lo:hi:
        %   R.mean_width     mean run width in pixels
        %   R.mean_runs      mean number of runs per row
        %   R.frac_two_plus  fraction of ALL rows with two or more runs
        %   R.rows_with_run  number of rows with at least one run
            [nR, ~] = size(rom_edge);
            widths = [];
            nruns  = zeros(nR, 1);
            for y = 1:nR
                cols = find(rom_edge(y, lo:hi));
                if isempty(cols), continue; end
                cols = cols(:)';
                brk  = [0, find(diff(cols) > 2), numel(cols)];
                for c = 1:numel(brk) - 1
                    seg = cols(brk(c)+1 : brk(c+1));
                    widths(end+1) = seg(end) - seg(1) + 1; %#ok<AGROW>
                end
                nruns(y) = numel(brk) - 1;
            end
            has = nruns > 0;
            if isempty(widths)
                R.mean_width = NaN;
            else
                R.mean_width = mean(widths);
            end
            if any(has)
                R.mean_runs = mean(nruns(has));
            else
                R.mean_runs = NaN;
            end
            R.frac_two_plus = 100 * nnz(nruns >= 2) / nR;
            R.rows_with_run = nnz(has);
        end

    % ====================================================================
    %  CHAPTER 6  Regression gate and file input and output
    %  compare_to_baseline treats a missing baseline as a failure, never as
    %  a skip, so a check that cannot run cannot pass.
    % ====================================================================

        function ok = compare_to_baseline(baseline_files, header, rows, keycols, numcols, tol)
        % COMPARE_TO_BASELINE  Regression check against a previously saved CSV.
        %
        %   baseline_files  cell array of candidate paths; the first that exists is used
        %   header, rows    new results in write_csv format
        %   keycols         columns identifying a row, e.g. {'Image','Variant'}
        %   numcols         columns to compare numerically
        %   tol             absolute tolerance
        %
        % Returns true only if a baseline exists AND every baseline row is
        % reproduced. A missing baseline is a FAILURE, not a skip: a check that
        % silently passes when it has nothing to compare against gives no protection.
            import rosed.*
            ok = true;
            file = '';
            for i = 1:numel(baseline_files)
                if exist(baseline_files{i}, 'file') == 2
                    file = baseline_files{i};
                    break;
                end
            end
            if isempty(file)
                ok = false;
                fprintf('  [regression] FAIL: no baseline found. Looked for:\n');
                for i = 1:numel(baseline_files)
                    fprintf('                 %s\n', baseline_files{i});
                end
                fprintf('  Put the v1 CSV in ..\\results\\tables\\ and rerun.\n');
                return;
            end
            [bh, br] = read_csv_simple(file);
            kb = cellfun(@(c) find(strcmp(bh, c), 1), keycols);
            nb = cellfun(@(c) find(strcmp(bh, c), 1), numcols);
            kn = cellfun(@(c) find(strcmp(header, c), 1), keycols);
            nn = cellfun(@(c) find(strcmp(header, c), 1), numcols);

            keystr = @(v) iff_num2str(v);
            newkeys = cell(size(rows, 1), 1);
            for r = 1:size(rows, 1)
                newkeys{r} = strjoin(cellfun(keystr, rows(r, kn), 'UniformOutput', false), '|');
            end

            nmatch = 0; nfail = 0; nmiss = 0;
            for r = 1:size(br, 1)
                key = strjoin(cellfun(keystr, br(r, kb), 'UniformOutput', false), '|');
                j = find(strcmp(newkeys, key), 1);
                if isempty(j)
                    nmiss = nmiss + 1;
                    fprintf('  [regression] MISSING  %s\n', key);
                    continue;
                end
                rowok = true;
                for c = 1:numel(numcols)
                    a = str2double(br{r, nb(c)});
                    b = rows{j, nn(c)};
                    if ischar(b), b = str2double(b); end
                    if abs(a - b) > tol
                        rowok = false;
                        fprintf('  [regression] DIFF  %-24s %-16s old %.4f  new %.4f\n', ...
                                key, numcols{c}, a, b);
                    end
                end
                if rowok, nmatch = nmatch + 1; else, nfail = nfail + 1; end
            end
            ok = (nfail == 0) && (nmiss == 0);
            fprintf('  [regression] %s: %d match, %d differ, %d missing  (%s)\n', ...
                    ternary(ok, 'PASS', 'FAIL'), nmatch, nfail, nmiss, file);
        end

        function s = iff_num2str(v)
            if ischar(v)
                x = str2double(v);
                if ~isnan(x), s = sprintf('%g', x); else, s = strtrim(v); end
            else
                s = sprintf('%g', v);
            end
        end

        function out = ternary(c, a, b)
            if c, out = a; else, out = b; end
        end

        function [header, rows] = read_csv_simple(file)
        % READ_CSV_SIMPLE  Portable CSV reader. Returns the header and every row as
        % a cell array of strings. Handles files written by MATLAB writetable.
            fid = fopen(file, 'r');
            if fid < 0, error('Cannot open %s', file); end
            header = strtrim(strsplit(fgetl(fid), ','));
            rows = {};
            while true
                ln = fgetl(fid);
                if ~ischar(ln), break; end
                if isempty(strtrim(ln)), continue; end
                rows(end+1, :) = strtrim(strsplit(ln, ',')); %#ok<AGROW>
            end
            fclose(fid);
        end

        function write_csv(file, header, rows)
        % WRITE_CSV  Portable CSV writer (no table type required).
        %   header : cell array of column names
        %   rows   : cell array, one row per line; numbers or strings
            d = fileparts(file);
            if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
            fid = fopen(file, 'w');
            if fid < 0, error('Cannot open %s for writing', file); end
            fprintf(fid, '%s\n', strjoin(header, ','));
            for r = 1:size(rows, 1)
                parts = cell(1, size(rows, 2));
                for c = 1:size(rows, 2)
                    v = rows{r, c};
                    if ischar(v)
                        parts{c} = v;
                    elseif abs(v - round(v)) < 1e-12
                        parts{c} = sprintf('%d', round(v));
                    else
                        parts{c} = sprintf('%.6f', v);
                    end
                end
                fprintf(fid, '%s\n', strjoin(parts, ','));
            end
            fclose(fid);
        end

    % ====================================================================
    %  CHAPTER 7  Small utilities
    %  sar_stretch is for display only and must never be fed to a detector.
    % ====================================================================

        function g = norm01(X)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            mn = min(X(:));  mx = max(X(:));
            if mx > mn
                g = (X - mn) / (mx - mn);
            else
                g = zeros(size(X));
            end
        end

        function G = sar_stretch(A, mode)
        % SAR_STRETCH  Display stretch. DISPLAY ONLY — never feed into a detector.
        %
        %   sar_stretch(A)            log transform + 2nd/98th percentile stretch.
        %                             For linear-scale SAR (simulated intensity).
        %   sar_stretch(A, 'linear')  1st/99th percentile stretch, no log.
        %                             For UMBRA GEC uint8, which the provider has
        %                             already scaled for display; a log on top of
        %                             that makes the image look nearly binary.
            if nargin < 2, mode = 'log'; end
            A = double(A);
            if strcmp(mode, 'log')
                L = log10(A + eps);  plo = 0.02;  phi = 0.98;
            else
                L = A;               plo = 0.01;  phi = 0.99;
            end
            f = sort(L(isfinite(L)));
            lo = f(max(1, round(plo * numel(f))));
            hi = f(min(numel(f), round(phi * numel(f))));
            G = (L - lo) / max(hi - lo, eps);
            G = min(max(G, 0), 1);
            G(~isfinite(G)) = 0;
        end

        function f = new_figure(name)
        % NEW_FIGURE  Visible figure in MATLAB, hidden when running headless.
            if exist('OCTAVE_VERSION', 'builtin') || ~usejava('desktop')
                f = figure('Name', name, 'Color', 'w', 'Visible', 'off');
            else
                f = figure('Name', name, 'Color', 'w');
            end
        end

        function save_png(fig, file)
        % SAVE_PNG  Save a figure as PNG at 150 dpi. Works in MATLAB and Octave.
            d = fileparts(file);
            if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
            print(fig, file, '-dpng', '-r150');
        end

        function e = sobel_edge(I)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
            gx = conv2(I, [-1 0 1; -2 0 2; -1 0 1], 'same');
            gy = conv2(I, [-1 -2 -1; 0 0 0; 1 2 1], 'same');
            mag = gx.^2 + gy.^2;
            e = mag > 4 * mean(mag(:));
        end

        function out = dilate_disk(BW, r)
        % [shared library] Extracted unchanged from rosed_pipeline.m so every experiment
        % calls identical code. Do not edit without rerunning the regression check.
        % Replacement for imdilate(BW, strel('disk', r)).
            [dx, dy] = meshgrid(-r:r, -r:r);
            mask = (dx.^2 + dy.^2) <= r^2;
            [rr, cc] = find(mask);
            rr = rr - (r+1);
            cc = cc - (r+1);
            [R, C] = size(BW);
            out = false(R, C);
            for i = 1:numel(rr)
                sr = rr(i);  sc = cc(i);
                rs = max(1, 1-sr) : min(R, R-sr);
                cs = max(1, 1-sc) : min(C, C-sc);
                out(rs, cs) = out(rs, cs) | BW(rs+sr, cs+sc);
            end
        end

    end

end
