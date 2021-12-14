function [fig, handles]= paintClusters(clustered, template, numLevels, isRecursive)
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
for i =1:numLevels
    angle = 2*pi/2^i;
    polarPoints{i} = angle * (1:2^i) - angle / 2;
end


if isRecursive
    [handles, ~] = paintRecursive(pax, xl, yl, axpos, ydir, clustered, template, polarPoints, [0,0], 1 / (2 * numLevels), 0);
else
    handles = paintIterative(pax, xl, yl, axpos, ydir, clustered, template, polarPoints,  1 / (2 * numLevels));
end
pax.set('XLim',[-0.5, 0.5]);
pax.set('YLim',[-0.5, 0.5]);
end

function ret=paintIterative(pax, xl, yl, axpos, ydir, clustered, template,polarPoints, offsetR)
cnt = 1;
for c=1:size(clustered,1)
    cluster = clustered(c, :);
    level = c-1;
    if level == 0
        center = [[0,0]];
        ret{cnt}.handle = paintSingle(find( cluster == 1), template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + center), offsetR/4);
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
            ret{cnt}.handle = paintSingle(find(selection), template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + center), offsetR/4);
            ret{cnt}.level = level;
            ret{cnt}.id = k;
            ret{cnt}.parent = 2 ^ (level-1) + ceil(k/2) - 1;
            cnt = cnt + 1;
        end
    end
    
end
end

function [ret, polarPoints] = paintRecursive(pax, xl, yl, axpos, ydir, clustered, template,polarPoints, center,  offsetR, level)
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

ret.handle= paintSingle(clustered.rootIndices, template, convertToFigureSpace(axpos, ydir, xl, yl, [0.5,0.5] + newCenter), offsetR/4);
if numel(clustered.parts)>0
    [ret.parts{1}, polarPoints]  = paintRecursive(pax, xl, yl, axpos, ydir, clustered.parts{1}, template,  polarPoints, newCenter,  offsetR, level+1);
    [ret.parts{2}, polarPoints]  = paintRecursive(pax, xl, yl, axpos, ydir, clustered.parts{2}, template,  polarPoints, newCenter,  offsetR, level+1);
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

function ax =paintSingle(rootIndices, template,  center, axisSize)
ax(1) = axes('pos',[center(1) - axisSize,  center(2) - axisSize/2, axisSize, axisSize]);
ax(2) = axes('pos',[center(1) ,  center(2) - axisSize/2, axisSize, axisSize]);
view(ax(1), 90, 0);
view(ax(2), -90, 0);
for i=1:2
    ax(i).Visible = false;
    ax(i).Color= 'none';
    shape =  shape3D;
    shape.Vertices = template.Vertices;
    shape.VertexRGB = 0.5 * ones(size(shape.Vertices));
    shape.VertexRGB(rootIndices, 2) = 1;
    shape.Faces = template.Faces;
    shape.Material = 'Dull';
    shape.ViewMode = 'solid';
    shape.ColorMode = 'texture';
    shape.RenderAxes = ax(i);
    shape.Visible = true;
    light = camlight(ax(i),'headlight');
    set(light,'Position',get(ax(i),'CameraPosition'));
end
end
