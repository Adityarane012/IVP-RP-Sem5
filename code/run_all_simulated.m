function run_all_simulated()
% RUN_ALL_SIMULATED  Runs every simulated experiment in order, logging each.
%
%   cd code
%   run_all_simulated
%
% Stops immediately if the regression check in rosed_pipeline.m or
% fig6e_extended.m fails, because every later result depends on them.
% Logs go to ..\results\logs\<script>.txt.
%
% Scripts run in the base workspace (evalin), so their "clear" commands
% cannot wipe this runner's loop.

    steps = { 'rosed_pipeline',        true;    % regression-gated
              'fig6e_extended',        true;    % regression-gated
              'fig6e_robustness',      false;
              'scale_sets_experiment', false;
              'fig6c_guardcells',      false;
              'classical_baselines',   false;
              'morphology_variant',    false };

    logdir = fullfile('..', 'results', 'logs');
    if ~exist(logdir, 'dir'), mkdir(logdir); end
    t_all = tic;

    for i = 1:size(steps, 1)
        name = steps{i, 1};
        logf = fullfile(logdir, [name '.txt']);
        if exist(logf, 'file'), delete(logf); end
        fprintf('\n##### [%d/%d] %s #####\n', i, size(steps, 1), name);
        t = tic;
        diary(logf);
        try
            evalin('base', 'REGRESSION_OK = true;');
            evalin('base', name);
        catch err
            diary off;
            fprintf(2, '\n%s FAILED: %s\nSee %s\n', name, err.message, logf);
            return;
        end
        diary off;
        fprintf('##### %s done in %.0f s #####\n', name, toc(t));

        if steps{i, 2} && ~evalin('base', 'REGRESSION_OK')
            fprintf(2, '\nSTOPPED: regression check failed in %s.\n', name);
            fprintf(2, 'Earlier results no longer reproduce. Do not continue until\n');
            fprintf(2, 'this is understood. Compare *_NEW.csv with the originals.\n');
            return;
        end
        close all;
    end
    fprintf('\nAll %d simulated experiments finished in %.1f min.\n', ...
            size(steps, 1), toc(t_all) / 60);
    fprintf('Tables: ..\\results\\tables   Figures: ..\\results\\figures   Logs: ..\\results\\logs\n');
end
