function [fig, handles]= paintClusters(clustered, template, numLevels, saveDir, isRecursive)
if nargin < 4
    isRecursive=true;
end
fig=figure();
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

levelShapes = cell(numLevels, 1);
if isRecursive
    [handles, ~, levelShapes] = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered, template, polarPoints, [0,0], 1 / (2 * numLevels), 0, levelShapes);
else
   [handles , levelShapes]= paintIterative(fig, pax, xl, yl, axpos, ydir, clustered, template, polarPoints,  1 / (2 * numLevels), levelShapes);
end
pax.set('XLim',[-0.5, 0.5]);
pax.set('YLim',[-0.5, 0.5]);
hold off;
f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.1, 0.8*(1/numLevels), 0.8];
t = tiledlayout(numLevels,2);
t.TileSpacing = 'none';
t.Padding = 'tight';
for i =1:numLevels
    shape = levelShapes{i};    
    ax = showPainted(f, shape, nan, nan, [nexttile(t), nexttile(t)]);
    daspect(ax(1), [1 1 1]);
    daspect(ax(2), [1 1 1]);
end
print(f,'-dsvg','-r300',[saveDir '/level' num2str(i-1)]);
end

function [ret, levelShapes]=paintIterative(fig, pax, xl, yl, axpos, ydir, clustered, template,polarPoints, offsetR, levelShapes)
cnt = 1;
ret = cell(2 ^ size(clustered,1)-1,1);
for c=1:size(clustered,1)
    cluster = clustered(c, :);
    level = c-1;
    if level == 0
        center = [[0,0]];
        pInds = find( cluster == 1);
        ret{cnt}.handle = paintSingle(fig, pInds, template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + center), offsetR/4);
        ret{cnt}.level = level;
        ret{cnt}.id = 1;
        ret{cnt}.parent = 0;
        cnt = cnt + 1;
    else
        levelPolarPoints = polarPoints{level};
        for k = 1:length(levelPolarPoints)
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
            ret{cnt}.handle = paintSingle(fig, pInds, template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + center), offsetR/4);
            levelShapes{level} = paintSingle(nan, pInds, template, nan, nan, levelShapes{level});
            ret{cnt}.level = level;
            ret{cnt}.id = k;
            ret{cnt}.parent = 2 ^ (level-1) + ceil(k/2) - 1;
            cnt = cnt + 1;
        end
    end

end
end

function [ret, polarPoints, levelShapes] = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered, template,polarPoints, center,  offsetR, level, levelShapes)
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

ret.handle= paintSingle(fig, clustered.rootIndices, template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + newCenter), offsetR/4);
if level >= 1
    levelShapes{level} = paintSingle(nan, clustered.rootIndices, template, nan, nan, levelShapes{level});
end
if numel(clustered.parts)>0
    [ret.parts{1}, polarPoints, levelShapes]  = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered.parts{1}, template,  polarPoints, newCenter,  offsetR, level+1, levelShapes);
    [ret.parts{2}, polarPoints, levelShapes]  = paintRecursive(fig, pax, xl, yl, axpos, ydir, clustered.parts{2}, template,  polarPoints, newCenter,  offsetR, level+1, levelShapes);
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

function ax =paintSingle(fig, rootIndices, template,  center, axisSize, shape)
if nargin < 6 || ~isa(shape,'shape3D')
    shape =  shape3D;
    shape.Vertices = template.Vertices;
    shape.VertexValue = zeros(shape.nVertices, 1);
    shape.VertexValue(rootIndices) = 1;
    shape.Faces = template.Faces;
    shape.Material = 'Dull';
    shape.ViewMode = 'solid';
    shape.ColorMode = 'indexed';
else
    shape.VertexValue(rootIndices) = max(shape.VertexValue) + 1;
end
if isa(fig,'matlab.ui.Figure')
    ax = showPainted(fig, shape, center,axisSize);
else
    ax = shape;
end

end

function ax=showPainted(fig, shape, center, axisSize, ax)
if nargin < 5
    if isnan(center)
        ax(1) = subplot(1,2,1);
        ax(2) = subplot(1,2,2);
    else
        ax(1) = axes(fig, 'pos',[center(1) - axisSize,  center(2) - axisSize/2, axisSize, axisSize]);
        ax(2) = axes(fig, 'pos',[center(1) ,  center(2) - axisSize/2, axisSize, axisSize]);
    end
end
view(ax(1), 90, 0);
view(ax(2), -90, 0);
maxlims = max(shape.Vertices);
minlims = min(shape.Vertices);

for i=1:2
    ax(i).Visible = false;
    ax(i).Color= 'none';
    s = clone(shape);
    s.ColorMode = 'indexed';
    s.RenderAxes = ax(i);
    s.Visible = true;
    s.PatchHandle.FaceColor = 'flat';
    colormap(ax(i),'colorcube');
    light = camlight(ax(i),'headlight');
    set(light,'Position',get(ax(i),'CameraPosition'));
end
xlim(ax(1), [ 1.3 * (maxlims(1) + minlims(1)) / 3, maxlims(1)])
end
