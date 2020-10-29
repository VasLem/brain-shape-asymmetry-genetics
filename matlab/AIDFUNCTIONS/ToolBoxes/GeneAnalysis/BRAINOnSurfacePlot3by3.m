function out = BRAINOnSurfacePlot3by3(rend,SegmentVal,levels,arrange)
    % First convert segmentvalues to surface values
    [~,VertexVal] = SegmentVal2VertexVal(rend,SegmentVal);
    out.f = figure;
    out.f.Color = [1 1 1];
    L = length(levels);
    for i=1:L
        %i=1;
       out.ax{i} = subplot(arrange(1),arrange(2),i);
       values = VertexVal(3,:);
       peer = out.ax{i};
       
    end
        
    
    
    
    






























end