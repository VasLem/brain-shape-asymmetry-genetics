function out = reduceTriangles(obj,R,varargin)
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

    if nargout == 1
       obj = clone(obj);
       %obj.Visible = false;
       out = obj;
    end
    [Vindex,Findex] = getVindexFindex(obj,varargin{:});
    if isempty(Findex), return; end
    
    % call matlab routine
    reduced = reducepatch(obj.Faces(:,Findex)',obj.Vertices',R);
    % relate reduced vertex indexing with old one
    [C,IA,IB] = intersect(obj.Vertices',reduced.vertices,'rows');
    [TF,LOC] = ismember(reduced.faces,IB);
    index = find(LOC);
    reduced.faces(index) = IA(LOC(index));
    % remove old faces and add new
    obj.Faces = [obj.Faces(:,setdiff((1:size(obj.Faces,2)),Findex)) reduced.faces'];
    % remove old vertices (not in faces)
    Vindex = Findex2Vindex(obj,(1:size(obj.Faces,2)));
    crop(obj,'VertexIndex',Vindex);
    
 
end