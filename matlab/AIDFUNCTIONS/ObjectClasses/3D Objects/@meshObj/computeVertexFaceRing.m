function vfring = computeVertexFaceRing(obj)

% computeVertexFacering - compute the faces adjacent to each vertex
%
%   ring = compute_vertex_face_ring(face);
%
%   Copyright (c) 2007 Gabriel Peyr?


%     if isempty(varargin)
%        face = obj.Faces;
%     else
%        face = obj.Faces(:,varargin{1}); 
%     end
    face = obj.Faces;
    nfaces = size(face,2);
    nverts = max(face(:));
    %nverts = length(unique(face(:)));

    vfring{nverts} = [];

    for i=1:nfaces
        for k=1:3
            vfring{face(k,i)}(end+1) = i;
        end
    end
end