function out = quad2tri(obj)
        if nargout == 1
            obj = clone(obj);
            out = obj;
        end
        if ~size(obj.Faces,1) == 4, return; end
        Quad = obj.Faces;
        Tri1 = zeros(3,size(Quad,2));
        Tri2 = zeros(3,size(Quad,2));
        for k=1:1:size(Quad,2)
            Tri1(:,k) = Quad(1:3,k);
            Tri2(1,k) = Quad(1,k);
            Tri2(2:end,k) = Quad(3:end,k);
        end
        obj.Faces = [Tri1,Tri2];
end