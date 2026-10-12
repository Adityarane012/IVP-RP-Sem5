function show_panel(nr, nc, k, G, M, ttl)
% [paper_figures helper]
    ax = subplot(nr, nc, k);
    if isempty(M)
        rgb = cat(3, G, G, G);
    else
        r = G; g = G; b = G;
        r(M) = 1; g(M) = 0; b(M) = 0;
        rgb = cat(3, r, g, b);
    end
    image(rgb, 'Parent', ax);
    axis(ax, 'image');  axis(ax, 'off');
    style_title(ax, ttl);
end
