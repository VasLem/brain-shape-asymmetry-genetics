function varargout = transformRigid(obj,T)
         if ~size(T,1)==4||~size(T,2)==4, error('Incorrect Transformation Matrix'); end
         if nargout == 1
             obj = clone(obj);
             obj.Visible = false;
             varargout{1} = obj;
         end
         obj.Vertices = transformRigid(obj.Vertices,T);
         updateChildren(obj,'Transform Rigid',T);
end