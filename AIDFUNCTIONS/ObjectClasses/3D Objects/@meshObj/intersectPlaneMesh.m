function out = intersectPlaneMesh(obj,N,P)
%         een vlak door de oorsprong met als normale the Y-as (vlak
%         parallell aan het X, Z vlak
%          N = [0 1 0]';
%          P = [0 0 0]';
         % getting distances from points to plane
         Dist = (obj.Vertices(1,:)-repmat(P(1),1,obj.nrV))*N(1) +...
                (obj.Vertices(2,:)-repmat(P(2),1,obj.nrV))*N(2) +...
                (obj.Vertices(3,:)-repmat(P(3),1,obj.nrV))*N(3);
            
         signTri =  reshape(sign(Dist(obj.Faces(:))),3,obj.nrF);
         
         diff = find(~(signTri(1,:)==signTri(2,:)));
         P1 = obj.Vertices(:,obj.Faces(1,diff));
         D1 = Dist(:,obj.Faces(1,diff));
         P2 = obj.Vertices(:,obj.Faces(2,diff));
         D2 = Dist(:,obj.Faces(2,diff));
         out = (repmat(D2,3,1).*P1 - repmat(D1,3,1).*P2)./repmat(D2-D1,3,1);
         
         diff = find(~(signTri(1,:)==signTri(3,:)));
         P1 = obj.Vertices(:,obj.Faces(1,diff));
         D1 = Dist(:,obj.Faces(1,diff));
         P2 = obj.Vertices(:,obj.Faces(3,diff));
         D2 = Dist(:,obj.Faces(3,diff));
         out =  [out (repmat(D2,3,1).*P1 - repmat(D1,3,1).*P2)./repmat(D2-D1,3,1)];
         
         diff = find(~(signTri(2,:)==signTri(3,:)));
         P1 = obj.Vertices(:,obj.Faces(2,diff));
         D1 = Dist(:,obj.Faces(2,diff));
         P2 = obj.Vertices(:,obj.Faces(3,diff));
         D2 = Dist(:,obj.Faces(3,diff));
         out =  [out (repmat(D2,3,1).*P1 - repmat(D1,3,1).*P2)./repmat(D2-D1,3,1)];
         out = unique(out','rows')';
end