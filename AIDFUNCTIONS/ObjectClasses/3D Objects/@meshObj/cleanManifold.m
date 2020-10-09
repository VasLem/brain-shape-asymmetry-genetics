function out = cleanManifold(obj)
         if nargout == 1, obj = clone(obj);out = obj; end
         % remove nan vertices
         [i,j] = find(isnan(obj.Vertices));
         if ~isempty(j), crop(obj,'VertexIndex',unique(j),'Action','delete'); end
         % remove vertices not part in triangles
         crop(obj,'VertexIndex',unique(obj.Tri(:)));
end