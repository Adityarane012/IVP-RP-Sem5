%% test_ordfilt2.m
% Confirms that ordfilt2 gives exactly the same output as our
% ordstat_filter for every domain the pipeline uses. Run from code/.
% Only replace ordstat_filter if every line prints PASS.

clear; clc;
import rosed.*          % the shared library, see code\rosed.m
rng(42);
P = randn(256).^2 + randn(256).^2;      % exponential clutter, no dependencies

allpass = true;
for s = [8 14 20]
    h = s/2;
    domL = false(s,s);  domL(:, 1:h)     = true;   % ROM left half
    domR = false(s,s);  domR(:, h+1:end) = true;   % ROM right half
    domF = true(s,s);                              % OS-CFAR full window

    tests = { 'ROM left',  domL, round(s*h/2);
              'ROM right', domR, round(s*h/2);
              'CFAR full', domF, round(0.75*s*s) };

    for t = 1:3
        a = ordstat_filter(P, tests{t,3}, tests{t,2});
        b = ordfilt2(P, tests{t,3}, tests{t,2}, 'symmetric');
        d = max(abs(a(:) - b(:)));
        ok = d == 0;
        allpass = allpass && ok;
        fprintf('s=%2d  %-10s  order=%4d  maxdiff=%g  %s\n', ...
                s, tests{t,1}, tests{t,3}, d, string(ok).replace(["true","false"],["PASS","FAIL"]));
    end
end

if allpass
    fprintf('\nALL PASS. Safe to switch ordstat_filter to ordfilt2.\n');
else
    fprintf('\nFAIL. Do not switch. Post this output.\n');
end
