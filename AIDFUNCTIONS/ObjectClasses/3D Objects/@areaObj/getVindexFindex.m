function [Vindex,Findex] = getVindexFindex(obj,varargin)
% tests the input arguments of a function in varargin whether Face index or
% Vertex index are given
% using a given Vertex index
    Input = find(strcmp(varargin, 'VertexIndex'));
    if isempty(Input)
        Vindex = (1:1:size(obj.Vertices,2));
    else
        Vindex = varargin{Input+1};
    end
    Findex = Vindex2Findex(obj,Vindex);
% using a given face index (OVERIDES A GIVEN VERTEX INDEX!!!!!)
    Input = find(strcmp(varargin, 'FaceIndex'));
    if ~isempty(Input)
        Findex = varargin{Input+1};
        Vindex = Findex2Vindex(obj,Findex);
    end
end