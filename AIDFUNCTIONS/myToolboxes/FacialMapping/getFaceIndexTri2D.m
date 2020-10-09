function F = getFaceIndexTri2D(QP,tri)
         nQP = size(QP,1);
         nF = size(tri.ConnectivityList,1);
         F = zeros(nQP,1);  
         for f=1:nF
             face = tri.ConnectivityList(f,:);
             IN = inpolygon(QP(:,1),QP(:,2),tri.Points(face,1),tri.Points(face,2));
             F(IN) = f;
         end
end