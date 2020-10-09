function [areaindex,p] = select3DRadiusArea(obj)
% selects a connected in terms of triangles area
areaindex = [];
%setptr(obj.Figure,'watch');
[p, v, vi, f, fi]  = select3DPoint(obj);
if isempty(p), return; end
distances = sqrt(sum((repmat(p,1,size(obj.CurrentMesh.Location,2))-obj.CurrentMesh.Location).^2));
areaindex = find(distances<=obj.SelectionSphere.Radius);
%setptr(obj.Figure,'arrow');
end
