function out = BRAINOnSurfaceArrange(rend,SegmentVal,levels,arrange,map,clim)
    % First convert segmentvalues to surface values
    [~,VertexVal] = SegmentVal2VertexVal(rend,SegmentVal);
    out.f1 = figure;
    out.f1.Color = [1 1 1];
    L = length(levels);
    for i=1:L
        %i=1;
       out.ax1{i} = subplot(arrange(1),arrange(2),i);
       renderBrainSurface(rend,VertexVal(levels(i),:),out.ax1{i});
       view(out.ax1{i},rend.viewval(1),0);
       light = camlight(out.ax1{i},'headlight');
       set(light,'Position',get(out.ax1{i},'CameraPosition'));
       drawnow;
       colormap(out.ax1{i},map);
       set(out.ax1{i},'clim',clim);
    end
    out.f2 = figure;
    out.f2.Color = [1 1 1];
    L = length(levels);
    for i=1:L
        %i=1;
       out.ax2{i} = subplot(arrange(1),arrange(2),i);
       renderBrainSurface(rend,VertexVal(levels(i),:),out.ax2{i});
       view(out.ax2{i},-1*rend.viewval(1),0);
       light = camlight(out.ax2{i},'headlight');
       set(light,'Position',get(out.ax2{i},'CameraPosition'));
       drawnow;
       colormap(out.ax2{i},map);
       set(out.ax2{i},'clim',clim);
    end
end