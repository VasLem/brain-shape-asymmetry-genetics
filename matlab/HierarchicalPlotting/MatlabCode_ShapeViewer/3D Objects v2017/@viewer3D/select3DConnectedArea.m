function areaindex = select3DConnectedArea(obj)
    % selects a connected in terms of triangles area
    areaindex = [];
    %[p,~,~,~,fi]  = select3DPoint(obj);
    p = mySelect3DPoint(obj);
    %if isempty(p), mySelect3DPoint(obj,'areaselection');end 
    if isempty(p), return; end
    vi = knnsearch(obj.CurrentShape.Vertices,p','K',1);
    distances = intraDistances(obj.CurrentShape,'VertexIndex',vi);
    areaindex = find(~isnan(distances));
    obj.SelectionDistances = distances;
    obj.SelectionTH = max(distances);
end
