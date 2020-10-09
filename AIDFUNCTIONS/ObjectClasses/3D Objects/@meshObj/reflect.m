function out = reflect(obj,Direction)
            if nargout == 1, obj = clone(obj);out = obj; end
            obj.Vertices(Direction,:) = -1*obj.Vertices(Direction,:);
            faces = obj.Faces;
            obj.Faces(2,:) = faces(3,:);
            obj.Faces(3,:) = faces(2,:);
            if ~isempty(obj.PoseLM)      
                obj.PoseLM.Vertices(Direction,:) = -1*obj.PoseLM.Vertices(Direction,:);
                if obj.PoseLM.nrV == 5
                   disp('Mirror 5');
                   obj.PoseLM.Vertices = obj.PoseLM.Vertices(:,[2 1 3 5 4]); 
                elseif obj.PoseLM.nrV == 7
                   disp('Mirror 7'); 
                   obj.PoseLM.Vertices = obj.PoseLM.Vertices(:,[2 1 3 5 4 7 6]);
                elseif obj.PoseLM.nrV == 12
                   disp('Mirror 12'); 
                   obj.PoseLM.Vertices = obj.PoseLM.Vertices(:,[4 3 2 1 7 6 5 10 9 8 11 12]);
                elseif obj.PoseLM.nrV == 14
                   disp('Mirror 14'); 
                   obj.PoseLM.Vertices = obj.PoseLM.Vertices(:,[4 3 2 1 7 6 5 10 9 8 11 12 14 13]);   
                else
                    orig = clone(obj.PoseLM);        
                    [N,D] = knn(kde(obj.PoseLM.Vertices,5),orig.Vertices,1);
                    obj.PoseLM.Vertices = obj.PoseLM.Vertices(:,N);
                end            
            end     
end