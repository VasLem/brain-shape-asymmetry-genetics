function varargout = transformRBF(obj,T)
         if nargout == 1
             obj = clone(obj);
             obj.Visible = false;
             varargout{1} = obj;
         end
         obj.Vertices = transformRBF(obj.Vertices,T);
         updateChildren(obj,'Transform RBF',T);
end