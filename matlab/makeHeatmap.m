
function fig=makeHeatmap(ret,traits,isPvalue, figPos)
    clrLim = [nanmin(abs(ret),[],'all'),nanmax(abs(ret),[],'all')];
    if nargin < 4
        figPos = nan;
    end
    fig=figure(Position=[0,0,400,1000]);
    [nr,nc] = size(ret);
    xrange = [1 nc]; % imagesc only needs the endpoints
    yrange = [1 nr];
    dx = diff(xrange)/(nc-1);
    dy = diff(yrange)/(nr-1);
    xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,nc+1);
    yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,nr+1);
    imagesc(xrange,yrange,abs(ret),'AlphaData',~isnan(ret)); hold on
    hm = mesh(xg,yg,zeros([nr,nc]+1));
    hm.FaceColor = 'none';
    hm.EdgeColor = 'k';
    ax = gca;
    t = num2cell(ret); % extact values into cells
    if isPvalue
        t = cellfun(@(x)sprintf('%.0e',x), t, 'UniformOutput', false); % convert to string
    else
        t = cellfun(@(x)sprintf('%.2f',x), t, 'UniformOutput', false); % convert to string
    end
    [y,x] = meshgrid(1:nr, 1:nc);
    % Draw Image and Label Pixels
    text(x(:), y(:), t', 'HorizontalAlignment', 'Center')
    
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
    a = colorbar();
    if isPvalue
        a.Label.String = 'P-value';
    else
        a.Label.String = 'Absol. Corr';
    end

    caxis(clrLim);
    set(ax,'xtick',1:length(traits));
    set(ax, 'YMinorTick','on')
%     set(ax,'YDir','normal')
    xticklabels(traits);
    if isPvalue
        set(ax,'ColorScale','log')
    end
    axis tight
    alignFigToContent(fig, ax, figPos);
end

function alignFigToContent(fig, h, figPos)
set(h, 'Units', 'pixels');
if isnan(figPos)
    figPos = get(h,'OuterPosition');
end
set(fig, 'Position', figPos);
set(h, 'Units',' normalized');
set(h, 'OuterPosition', [0,0,1,1]);
end