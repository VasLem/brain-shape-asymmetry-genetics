function [ax,pRet]=showPaintedDoubleFace(fig, shape, center, axisSize, ax, jetcmap, clims, invcolor,pToUse)
% Show supplied shape on the given figure, using a certain center and axis size if supplied
% on new axes or a set of two if supplied in the argumnet ax. Returns the
% modified axes. Adds jet cmap if jetcmap is true, setting colormap to jet instead of colorcube,
% and using the supplied clims. invcolor is the invalid color to be added
% to the hole, defaults to black if not given

%
if nargin < 5 || ~isa(ax,'matlab.graphics.axis.Axes')
    if isnan(center)
        %         ax(1) = subplot(1,2,1);
        %         ax(2) = subplot(1,2,2);
        ax = subplot(1,1,1);
    else
        ax = axes(fig, 'pos',[center(1) - axisSize/2,  center(2) - axisSize/2, axisSize, axisSize]);
        %         ax(1) = axes(fig, 'pos',[center(1) - axisSize,  center(2) - axisSize/2, axisSize, axisSize]);
        %         ax(2) = axes(fig, 'pos',[center(1) ,  center(2) - axisSize/2, axisSize, axisSize]);
    end
end
if nargin < 6 || isnan(jetcmap)
    jetcmap = false;
end
if nargin < 7 || any(isnan(clims))
    clims = [0,1];
end
if nargin < 8 | any(isnan(invcolor))
    invcolor = 'white';
end
if isnan(invcolor)
    invcolor = 'white';
end
if nargin < 9
    pToUse =nan;
end
ax.Visible=false;

for i=1:2
    s = clone(shape);
    if i==2
        s.Vertices(:,1) = - s.Vertices(:,1);
        s.Vertices(:, 2) = - s.Vertices(:, 2); % to show the left hemisphere
    end
    
    if i == 1
        s.Vertices(:,2) = s.Vertices(:,2) - abs(max(s.Vertices(:,2))) - 0.1;
    else
        s.Vertices(:,2) = s.Vertices(:,2) + abs(min(s.Vertices(:,2))) + 0.1;
    end

    s.Vertices(:,1) = s.Vertices(:,1) - mean(s.Vertices(:,1));
    s.ColorMode = 'indexed';
    if i == 1
           [px,py,pz] = addRectangle(ax,s,invcolor, pToUse);
           pRet = [px py pz];
           view(ax, 90, 0);
    end
    s.RenderAxes = ax;
    s.Visible = true;
    s.PatchHandle.FaceColor = 'flat';
    
    if i == 2
        light = camlight(ax,'headlight');
        set(light,'Position',get(ax,'CameraPosition'));
    end
end
daspect(ax, [1 1 1]);
axis(ax,'image');
axis(ax,'off');
end
function [x,y,z]=addRectangle(ax, shape, color, pToUse)
if isnan(pToUse)
    maxlims = max(shape.Vertices);
    x = maxlims(1)/30;
    y = double(shape.Vertices( abs(shape.Vertices(:,1) - x) < 0.001 * max(shape.Vertices(:,1)) ,2));
    z = double(shape.Vertices( abs(shape.Vertices(:,1) - x) < 0.001 * max(shape.Vertices(:,1)),3));
    p = boundary(y, z);
    y = y(p);
    z = z(p);
    x = repmat(x, length(p),1);
else
    x = pToUse(:,1);
    y = pToUse(:,2);
    z = pToUse(:,3);
end
patch(ax, x,y,z,color,'EdgeColor','none', 'FaceLighting','none');
end