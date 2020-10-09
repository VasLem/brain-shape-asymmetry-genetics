function Vindex = Findex2Vindex(obj,Findex)
% converts a Face Index (selection of mesh faces) into Vertex index
% (selection of mesh vertices)
    faces = obj.Faces(:,Findex);
    Vindex = unique(faces(:));
end