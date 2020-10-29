function Findex = Vindex2Findex(obj,Vindex)
% converts a Vertex Index (selection of mesh vertices) into Face index
% (selection of mesh faces)
    tmp = ismember(obj.Faces,Vindex);
    [i,j] = find(tmp==0);
    Findex = setdiff((1:size(obj.Faces,2)),j);

end