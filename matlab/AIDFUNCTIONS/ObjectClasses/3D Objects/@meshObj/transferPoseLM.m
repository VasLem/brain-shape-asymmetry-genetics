function out = transferPoseLM(obj,refobj)
         if nargout == 1
            obj = clone(obj);
            out = obj;
         end
         delete(obj.PoseLM);
         bar = cart2bary(refobj.Vertices,refobj.Faces,refobj.PoseLM.Vertices,refobj.PoseLM.Fi);
         cart = bary2cart(obj.Vertices,obj.Faces,bar,refobj.PoseLM.Fi);
         obj.PoseLM = LMObj('Vertices',cart,'Fi',refobj.PoseLM.Fi);
end