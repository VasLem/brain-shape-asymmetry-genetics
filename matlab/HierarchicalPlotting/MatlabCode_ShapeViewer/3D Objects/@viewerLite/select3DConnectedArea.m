function areaindex = select3DConnectedArea(obj)
% selects a connected in terms of triangles area
areaindex = [];
[p, v, vi, f, fi]  = select3DPoint(obj);
if isempty(p), return; end
distances = intraDistances(obj.CurrentMesh,'VertexIndex',obj.CurrentMesh.Tri(:,fi)');
areaindex = find(~isnan(distances));
obj.UserData.Distances = distances;
obj.UserData.TH = max(distances);
end
