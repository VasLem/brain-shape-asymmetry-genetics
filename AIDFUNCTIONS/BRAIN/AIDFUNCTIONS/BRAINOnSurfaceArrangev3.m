function out = BRAINOnSurfaceArrangev3(in,rend,SegmentVal,levels,map,clim)
    % First convert segmentvalues to surface values
    [~,VertexVal] = SegmentVal2VertexVal(rend,SegmentVal);
    out = in;
    L = length(levels);
    counter = 0;
    titlenames = {'i' 'ii' 'iii' 'iv' 'v' 'vi' 'vii' 'viii' 'ix'};
    for i=1:L
        %i=1;
       counter = counter+1;
       %out.ax1{i} = subplot(arrange(1),arrange(2),counter);
       colormap(out.ax1{i},map);
       if ~isempty(clim), set(out.ax1{i},'clim',clim);end
       renderBrainSurface(rend,VertexVal(levels(i),:),out.ax1{i});
       view(out.ax1{i},rend.viewval(1),0);
       light = camlight(out.ax1{i},'headlight');
       set(light,'Position',get(out.ax1{i},'CameraPosition'));
       drawnow;
       title(out.ax1{i},titlenames(levels(i)));
       counter = counter+1;
       %out.ax2{i} = subplot(arrange(1),arrange(2),counter);
       renderBrainSurface(rend,VertexVal(levels(i),:),out.ax2{i});
       view(out.ax2{i},-1*rend.viewval(1),0);
       light = camlight(out.ax2{i},'headlight');
       set(light,'Position',get(out.ax2{i},'CameraPosition'));
       drawnow;
       colormap(out.ax2{i},map);
       if ~isempty(clim), set(out.ax2{i},'clim',clim);end
    end
end