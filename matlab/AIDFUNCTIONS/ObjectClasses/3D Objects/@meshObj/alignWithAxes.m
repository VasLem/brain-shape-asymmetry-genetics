function obj = alignWithAxes(obj)
         if nargout==1, obj = clone(obj); end
         % Step 1, align centroid with axes origin
            C = centroid(obj);obj.Vertices = obj.Vertices-repmat(C,1,obj.nrV);
         % Step 2, align PCA axes with axes
            EigVec = pca(obj.Vertices);
            obj.Vertices = inv(EigVec)*obj.Vertices; %#ok<*MINV>
            C = centroid(obj);obj.Vertices = obj.Vertices-repmat(C,1,obj.nrV);
         % Step 3, make sure Y is height, X is width and Z is depth (FOR FACES)    
            R = inv([0 1 0;1 0 0;0 0 1]);obj.Vertices = R*obj.Vertices;
            C = centroid(obj);obj.Vertices = obj.Vertices-repmat(C,1,obj.nrV);   
end