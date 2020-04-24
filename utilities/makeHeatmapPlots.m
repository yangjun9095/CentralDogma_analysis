%%%%%%% Nested function to make basic and PBoC-style heatmap plots %%%%%%%
function heatmap_fig = makeHeatmapPlots(Image, VisibleOn, Title, CLabel,Colormap_heat,PixelSize,lb,ub)
    if VisibleOn
        heatmap_fig = figure;
    else
        heatmap_fig = figure('Visible','Off');
    end
    colormap(Colormap_heat)
    snippet_im = imagesc(Image);
    title(Title)
    caxis([lb ub])
    c = colorbar;
    c.Ticks = linspace(lb,ub,3);
    ylabel('position (\mum)','FontSize',15)
    xlabel('position (\mum)','FontSize',15)
    ylabel(c,CLabel,'FontSize',15)
    xy_lim = (size(Image,1) + 1) / 2 * PixelSize; % um, calculated with center of center pixel as 0 
    set(gca,'xtick',([-1.0 0 1.0] + xy_lim) ./ PixelSize, ...
            'xticklabel',[-1.0, 0, 1.0])
    set(gca,'ytick',([-1.0 0 1.0] + xy_lim) ./ PixelSize, ...
            'yticklabel',[-1.0, 0, 1.0])

    % make paper fig in PBoC style
    StandardFigurePBoC(snippet_im,gca);
end