function f = paper_fig(name, pos)
% [paper_figures helper]
% White figure with black text, whatever the MATLAB desktop theme is.
    import rosed.*          % the shared library, see code\rosed.m
    f = new_figure(name);
    try, theme(f, 'light'); catch, end        % MATLAB R2025a and later
    set(f, 'Color', 'w', 'Position', pos);
end
