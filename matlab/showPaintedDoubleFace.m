function ax=showPaintedDoubleFace(fig, shape, center, axisSize, ax, jetcmap, clims)
% Show supplied shape on the given figure, using a certain center and axis size if supplied
% on new axes or a set of two if supplied in the argumnet ax. Returns the
% modified axes. Adds jet cmap if jetcmap is true, setting colormap to jet instead of colorcube,
% and using the supplied clims.
%
if nargin < 5
    if isnan(center)
        ax(1) = subplot(1,2,1);
        ax(2) = subplot(1,2,2);
    else
        ax(1) = axes(fig, 'pos',[center(1) - axisSize,  center(2) - axisSize/2, axisSize, axisSize]);
        ax(2) = axes(fig, 'pos',[center(1) ,  center(2) - axisSize/2, axisSize, axisSize]);
    end
end
if nargin < 6
    jetcmap = false;
end
if nargin < 7
    clims = [0,1];
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
    light = camlight(ax(i),'headlight');
    set(light,'Position',get(ax(i),'CameraPosition'));    
    if numel(jetcmap)==1
        if ~jetcmap
            colormap(ax(i),'colorcube');
            if length(unique(s.VertexValue))==1
                caxis(ax(i), [0,1]);
            end
        else
            colormap(ax(i), 'jet');
            caxis(ax(i), clims);
        end
    else
        s.ColorMode = 'texture';
    end
end
xlim(ax(1), [ 1.3 * (maxlims(1) + minlims(1)) / 3, maxlims(1)])
ylim(ax(1), [min(s.Vertices(:,2)), max(s.Vertices(:,2))] + 0.1);
ylim(ax(2), [min(s.Vertices(:,2)), max(s.Vertices(:,2))] - 0.1);
daspect(ax(1), [0.5 1 1]);
daspect(ax(2), [0.5 1 1]);
end

