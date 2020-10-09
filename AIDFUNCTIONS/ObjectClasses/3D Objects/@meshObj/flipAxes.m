function obj = flipAxes(obj,ax)
         if nargout==1, obj = clone(obj); end
         % Step 1, align centroid with axes origin
            C = centroid(obj);obj.Vertices = obj.Vertices-repmat(C,1,obj.nrV);
         % Step 2, flip directions asked for   
            flips = ones(1,3);flips(ax) = -1;
            R = inv(repmat(flips,3,1).*eye(3,3));obj.Vertices = R*obj.Vertices;
            C = centroid(obj);obj.Vertices = obj.Vertices-repmat(C,1,obj.nrV);
end