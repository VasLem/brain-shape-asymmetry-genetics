function out = FACEOnSurfaceArrange(peer,rend,SegmentVal,levels,arrange,map,clim)
    % First convert segmentvalues to surface values
    [~,VertexVal] = SegmentVal2VertexVal(rend,SegmentVal);
    out.f1 = peer;
    %figure(peer);
    out.f1.Color = [1 1 1];
    L = length(levels);
    titlenames = {'i' 'ii' 'iii' 'iv' 'v' 'vi' 'vii' 'viii' 'ix'};
    for i=1:L
        %i=1;
       out.ax1{i} = subplot(arrange(1),arrange(2),i,'Parent',peer);
       out.scan{i} = renderFaceSurface(rend,VertexVal(i,:),out.ax1{i});
       %view(out.ax1{i},rend.viewval(1),0);
       light = camlight(out.ax1{i},'headlight');
       set(light,'Position',get(out.ax1{i},'CameraPosition'));
       drawnow;
       colormap(out.ax1{i},map);
       title(out.ax1{i},titlenames(i));
       if ~isempty(clim),set(out.ax1{i},'clim',clim);end
    end
end