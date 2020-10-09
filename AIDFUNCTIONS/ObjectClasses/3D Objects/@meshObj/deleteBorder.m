function out = deleteBorder(obj,times)
         if nargin ==1, times = 1; end
         if nargout == 1
             obj = clone(obj);
             out = obj;
         end
         if times>1, deleteBorder(obj,times-1);return;end
         border(obj);
         crop(obj,'VertexIndex',obj.Border.VerticesIndex,'Action','delete');
end