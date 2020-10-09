function varargout = reduceTriangles(obj,R)
% REDUCETRIANGLES(obj, R) reduces the number of faces in patch P while trying
%     to preserve the overall shape of the patch.  If R is less than or
%     equal to 1, R is interpreted as a fraction of the original faces; for
%     example, if R is 0.2, 20% of the faces will be kept. If R is greater
%     than 1, then R is the target number of faces. For example, if R were
%     400, then the number of faces would be reduced until there were 400
%     faces remaining. If the patch contains non-shared vertices, shared
%     vertices are computed before reduction. If the faces of the patch are
%     not triangles, the faces are triangulated before reduction. The faces
%     returned are always triangles.
% varargin: 'VertexIndex', 'Faceindex'

    old_vertices = obj.Vertices;
    
    if nargout == 1
        % was already cloned with the mesh cloning
        varargout{1} = reduceTriangles(obj.Parent,R,'VertexIndex',obj.VerticesIndex);
    else
       reduceTriangles(obj.Parent,R,'VertexIndex',obj.VerticesIndex); 
       [C,IA,IB] = intersect(obj.Parent.Vertices',old_vertices','rows');
       obj.VerticesIndex = IA';
    end
end