function subdivideTriangles(obj,mode,val)
         nrV = size(obj.Parent.Vertices,2);
         subdivideTriangles(obj.Parent,mode,val,'VertexIndex',obj.VerticesIndex); 
         newnrV = size(obj.Parent.Vertices,2);
         obj.VerticesIndex = [obj.VerticesIndex (nrV+1:1:newnrV)];
end