function ax=showPaintedDoubleFace(fig, shape, center, axisSize, ax, jetcmap, clims, invcolor)
% Show supplied shape on the given figure, using a certain center and axis size if supplied
% on new axes or a set of two if supplied in the argumnet ax. Returns the
% modified axes. Adds jet cmap if jetcmap is true, setting colormap to jet instead of colorcube,
% and using the supplied clims. invcolor is the invalid color to be added
% to the hole, defaults to black if not given

%
if nargin < 5
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
if nargin < 6
    jetcmap = false;
end
if nargin < 7
    clims = [0,1];
end
if nargin < 8
    invcolor = 'white';
end
ax.Visible=false;
if numel(jetcmap)==1
    if ~jetcmap
        invcolor = colorcube(1);
    end
    if jetcmap
        invcolor = jet(1);
    end
end

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
           addRectangle(ax,s,invcolor);
           view(ax, 90, 0);
    end
    s.RenderAxes = ax;
    s.Visible = true;
    s.PatchHandle.FaceColor = 'flat';
    
    if i == 2

        light = camlight(ax,'headlight');
        set(light,'Position',get(ax,'CameraPosition'));
        if numel(jetcmap)==1
            if ~jetcmap
                colormap(ax,'colorcube');
             
                if length(unique(s.VertexValue))==1
                    caxis(ax, [0,1]);
                end
            else
                colormap(ax, 'jet');
                caxis(ax, clims);
            end
        end
    end
    if numel(jetcmap) ~= 1
        s.ColorMode = 'texture';
        colormap(ax,jetcmap)
    end


end

daspect(ax, [1 1 1]);
end
function addRectangle(ax, shape, color)
maxlims = max(shape.Vertices);
x = maxlims(1)/30;
y = double(shape.Vertices( abs(shape.Vertices(:,1) - x) < 0.001 * max(shape.Vertices(:,1)) ,2));
z = double(shape.Vertices( abs(shape.Vertices(:,1) - x) < 0.001 * max(shape.Vertices(:,1)),3));

p = boundary(y, z);
y = y(p);
z = z(p);
x = repmat(x, length(p),1);
p = patch(ax, x,y,z,color,'EdgeColor','none', 'FaceLighting','none');
end