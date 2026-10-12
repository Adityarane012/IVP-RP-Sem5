function save_paper(f, file)
% [paper_figures helper]
% 300 dpi, white background.
    try
        exportgraphics(f, file, 'Resolution', 300, 'BackgroundColor', 'white');
    catch
        print(f, file, '-dpng', '-r300');
    end
    fprintf('saved %s\n', file);
end
