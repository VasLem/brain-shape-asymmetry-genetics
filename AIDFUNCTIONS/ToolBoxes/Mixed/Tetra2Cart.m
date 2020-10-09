function [p,Tetra] = Tetra2Cart(C)
         Vertices = [0 0 0;0 1 0;0.866 0.5 0; 0.2887 0.5 0.866]';
         X = Vertices(1,:);
         Y = Vertices(2,:);
         Z = Vertices(3,:);
         Tri = [1 2 3;1 2 4;2 3 4;3 1 4]';
         Tetra = meshObj('Vertices',Vertices,'Faces',Tri);
         Tetra.ViewMode = 'Wireframe';
         Tetra.SingleColor = [1 1 1];        
         T = [X(1)-X(4) X(2)-X(4) X(3)-X(4);...
             Y(1)-Y(4) Y(2)-Y(4) Y(3)-Y(4);...
             Z(1)-Z(4) Z(2)-Z(4) Z(3)-Z(4)];
         R = T*C';
         R = R+repmat(Vertices(:,4),1,size(C,1));
         p = LMObj('Vertices',R);
end