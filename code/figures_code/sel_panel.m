function sel_panel(k, sel, cmap, ttl)
% [paper_figures helper]
    ax = subplot(1, 2, k);
    imagesc(sel, 'Parent', ax);
    set(ax, 'CLim', [-0.5 3.5], 'XColor', 'k', 'YColor', 'k', 'FontSize', 10);
    colormap(ax, cmap);  axis(ax, 'image');
    cb = colorbar(ax);
    labels = {'none', '8', '14', '20'};
    try
        set(cb, 'Ticks', 0:3, 'TickLabels', labels, 'Color', 'k');    % MATLAB
    catch
        set(cb, 'ytick', 0:3, 'yticklabel', labels);                  % Octave
    end
    style_title(ax, ttl);
end
