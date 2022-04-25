function [fig, fig2, handles]= paintClusters(clustered, template, numLevels, isRecursive, values, invcolor,map, tilesLayout)
% Paint clusters with a certain number of levels numLevels of a shape3d
% template using concentric circles.
if nargin < 4
    isRecursive=false;
end
if nargin < 5
    values = [];
end
if nargin < 6
    invcolor = 'white';
end
if nargin < 7
    map = 'parula';
end

if ~isempty(values) && isRecursive
        error("Cannot visualize custom values when painting recursively.")
end
fig=figure();
fig.Units = 'normalized';
pax = gca();
hold on;
pax.Visible = false;
xl = xlim();
yl = ylim();
axpos = get(pax,'position');
ydir = get(pax, 'ydir');
polarPoints = cell(numLevels,1);
for i =1:numLevels
    angle = 2*pi/2^i;
    polarPoints{i} = angle * (1:2^i) - angle / 2;
end

levelShapes = cell(numLevels + 1, 1);
if isRecursive
    [handles, ~, levelShapes] = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered, template, polarPoints, [0,0], 1 / (2 * numLevels), 0, levelShapes,invcolor,map);
else
   [handles , levelShapes]= paintIterative(fig, pax, xl, yl, axpos, ydir, clustered, template, polarPoints,  1 / (2 * numLevels), levelShapes,values,invcolor,map);
end
pax.set('XLim',[-0.5, 0.5]);
pax.set('YLim',[-0.5, 0.5]);
hold off;
fig2 = figure;
fig2.Units = 'normalized';
fig2.Position = [0.1,0.1, 0.8*(1/numLevels), 0.8];
if nargin < 8 || ~strcmp(tilesLayout, 'square')
    if ~isempty(values)
        % give space to colorbar
        t = tiledlayout(numLevels + 2,1);
    else
        t = tiledlayout(numLevels + 1,1);
    end
    t.TileSpacing = 'none';
    t.Padding = 'tight';
else
    if strcmp(tilesLayout,'square')
        fig2.Position = [0.1,0.1, 0.6, 0.4];
        if ~isempty(values)
            t = tiledlayout(3,2);
        else
            t = tiledlayout(2,2);
        end
    end
end

minValue = 0;
maxValue = 1;
if ~isempty(values)
minValue = min(values);
maxValue = max(values);
end
for i =1:(numLevels + 1)
    if isempty(values)
        if i == 1
            continue
        end
    end
    shape = levelShapes{i};    
    ax = showPaintedDoubleFace(fig2, shape, nan, nan, nexttile(t), ~isempty(values),[minValue,maxValue], invcolor);
    if i == 1
        title(ax,'root','FontSIze',16)
    else
        title(ax,['level ' num2str(i-1)],'FontSize',16)
    end
    daspect(ax, [1 1 1]);
end

if ~isempty(values)
    h = axes(fig2,'visible','off');
    h.Units = 'normalized';
    h.Position = [0.1,0,0.8,0.7];
    ticks = linspace(minValue,maxValue, 10);
    colorbar(h, 'South','TickLabels',  round(ticks*1000)/1000, ...
        'Ticks',linspace(minValue,maxValue, 10),'FontWeight','bold');
    colormap('jet');
    caxis(h,[minValue, maxValue]);
end
set(fig2, 'InvertHardCopy', 'off');
set(fig2, 'Color', 'white');
end

function [ret, levelShapes]=paintIterative(fig, pax, xl, yl, axpos, ydir, clustered, template,polarPoints, offsetR, levelShapes, values,invcolor,map)
cnt = 1;
ret = cell(2 ^ size(clustered,1)-1,1);
for c=1:size(clustered,1)
    cluster = clustered(c, :);
    if ~isempty(values)
        value = values(cnt);
    else
        value = nan;
    end
    level = c-1;
    if level == 0
        center = [[0,0]];
        pInds = find( cluster == 1);
        ret{cnt}.handle = paintSingle(fig, pInds, template, ...
            convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + center), offsetR, nan, value,invcolor,map);

        levelShapes{level + 1} = paintSingle(nan, pInds, template, nan, nan, levelShapes{level + 1},value,invcolor,map);
        ret{cnt}.level = level;
        ret{cnt}.id = 1;
        ret{cnt}.parent = 0;
        cnt = cnt + 1;
    else
        levelPolarPoints = polarPoints{level};
        for k = 1:length(levelPolarPoints)
            if ~isempty(values)
                value = values(cnt);
            else
                value = nan;
            end
            angle = levelPolarPoints(k);
            if level == 1
                prevAngle = 0;
            else
                prevAngle = polarPoints{level-1}(ceil(k/2));
            end
            selection = cluster == k;
            if ~any(selection), cnt=cnt+1; continue; end
            [prevX, prevY] = pol2cart(prevAngle, (level-1)*offsetR);
            [x, y] = pol2cart(angle, level * offsetR);
            line(pax, [prevX, x], [prevY, y],'Color','black');
            set(pax,'xlim',[0,1]);
            set(pax,'ylim',[0,1]);
            center = [x,y];
            pInds = find(selection);
            ret{cnt}.handle = paintSingle(fig, pInds, template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + center), ...
                offsetR, nan, value,invcolor,map);
            levelShapes{level + 1} = paintSingle(nan, pInds, template, nan, nan, levelShapes{level + 1},value,invcolor,map);
            ret{cnt}.level = level;
            ret{cnt}.id = k;
            ret{cnt}.parent = 2 ^ (level-1) + ceil(k/2) - 1;
            cnt = cnt + 1;
        end
    end

end
end

function [ret, polarPoints, levelShapes] = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered, template,polarPoints, center,  offsetR, level, levelShapes, invcolor,map)
if level > 0
    levelPolarPoints = polarPoints{level};
    angle = levelPolarPoints(end);
    polarPoints{level} = levelPolarPoints(1:end-1);
    [x, y] = pol2cart(angle, level * offsetR);
    newCenter = [x,y];
    line(pax, [center(1), newCenter(1)], [center(2), newCenter(2)],'Color','black');
    set(pax,'xlim',[0,1]);
    set(pax,'ylim',[0,1]);
else
    newCenter = [0,0];
end

ret.handle= paintSingle(fig, clustered.rootIndices, template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + newCenter), ...
    offsetR, nan, nan,invcolor,map);
levelShapes{level + 1} = paintSingle(nan, clustered.rootIndices, template, nan, nan, levelShapes{level + 1},nan,invcolor,map);
if numel(clustered.parts)>0
    [ret.parts{1}, polarPoints, levelShapes]  = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered.parts{1}, template,  polarPoints, newCenter,  offsetR, level+1, levelShapes,invcolor,map);
    [ret.parts{2}, polarPoints, levelShapes]  = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered.parts{2}, template,  polarPoints, newCenter,  offsetR, level+1, levelShapes,invcolor,map);
else
    if numel(polarPoints)>=level+1
        if numel(polarPoints{level+1}) > 0
            polarPoints{level+1} = polarPoints{level+1}(1:end-2);
        end
    end
end
end


function xycFigNorm = convertToFigureSpace(axpos, ydir, xl, yl, xyc)
xycNorm = (xyc - [xl(1),yl(1)])./[range(xl),range(yl)]; %normalized to axis
if strcmpi(ydir, 'reverse')
    xycNorm(2) = 1 - xycNorm(2); % This is needed iff y axis direction is reversed!
end
xycFigNorm = axpos(1:2) + axpos(3:4).*xycNorm; %normalized to figure
end

function ax =paintSingle(fig, rootIndices, template,  center, axisSize, shape, value,invcolor,map)
if  ~isa(shape,'shape3D')
    shape =  shape3D;
    shape.Vertices = template.Vertices;
    shape.VertexValue = zeros(shape.nVertices, 1);
    shape.SingleColor = [0,0,0];
    if isnan(value)
        value = 1;
    end
    shape.VertexValue(rootIndices) = value;
    shape.Faces = template.Faces;
    shape.Material = 'Dull';
    shape.ViewMode = 'solid';
    shape.ColorMode = 'indexed';
else
    if isnan(value)
        value =  max(shape.VertexValue) + 1;
    end
    shape.VertexValue(rootIndices) = value;
end
if isa(fig,'matlab.ui.Figure')
    ax = showPaintedDoubleFace(fig, shape, center,axisSize,nan,false,[0,1], invcolor);
    colormap(ax,map)
else
    ax = shape;
end
end


