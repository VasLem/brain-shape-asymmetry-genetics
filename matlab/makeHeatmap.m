
function fig=makeHeatmap(ret,traits,isPvalue)
    clrLim = [nanmin(ret,[],'all'),nanmax(ret,[],'all')];

    fig=figure(Position=[0,0,400,1000]);
    [nr,nc] = size(ret);
    xrange = [1 nc]; % imagesc only needs the endpoints
    yrange = [1 nr];
    dx = diff(xrange)/(nc-1);
    dy = diff(yrange)/(nr-1);
    xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,nc+1);
    yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,nr+1);
    imagesc(xrange,yrange,ret,'AlphaData',~isnan(ret)); hold on
    hm = mesh(xg,yg,zeros([nr,nc]+1));
    hm.FaceColor = 'none';
    hm.EdgeColor = 'k';
    ax = gca;
    t = num2cell(ret); % extact values into cells
    if isPvalue
        t = cellfun(@(x)sprintf('%.1e',x), t, 'UniformOutput', false); % convert to string
    else
        t = cellfun(@(x)sprintf('%.2f',x), t, 'UniformOutput', false); % convert to string
    end
    [y,x] = meshgrid(1:nr, 1:nc);
    % Draw Image and Label Pixels
    text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
    
    l = linspace(0.3,1,100);
    s = ones(size(l));
    if isPvalue
        h = 0.6 + zeros(size(l));
    else
        h = zeros(size(l));
    end
    map = hsl2rgb([h;s;l]');
    if isPvalue
        map = flip(map,1);
    end
    colormap(ax,map);
    colorbar();
    caxis(clrLim);
    set(ax,'xtick',1:length(traits));
    set(ax, 'YMinorTick','on')
    set(ax,'YDir','normal')
    xticklabels(traits);
    if isPvalue
        set(ax,'ColorScale','log')
    end
    axis tight
    alignFigToContent(fig, ax);
end

function alignFigToContent(fig, h)
set(h, 'Units', 'pixels');
pos = get(h,'OuterPosition');
set(fig, 'Position',pos);
set(h, 'Units',' normalized');
set(h, 'OuterPosition', [0,0,1,1]);
end