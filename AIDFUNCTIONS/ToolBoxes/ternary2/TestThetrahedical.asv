%% Tetrahedrial
figure; hold on;
X = [0 0 1 0;...
     0 0 0.5 0];
Y = [0 1 0.5 0;0 1 0.5 0];
Z = [0 0 0 0;0 0 1 0];

h = fill3(,,,'w','linewidth',2);
h = fill3([],[],[],'w','linewidth',2);
%h = fill3([0 0 1 0],[0 1 0.5 0],[0 0 0 0],'w','linewidth',2);


Vertices = [0 0 0;0 1 0;0.866 0.5 0; 0.2887 0.5 0.866]';
X = Vertices(1,:);
Y = Vertices(2,:);
Z = Vertices(3,:);
Tri = [1 2 3;1 2 4;2 3 4;3 1 4]';
Tetra = meshObj('Vertices',Vertices,'Faces',Tri);

v = viewer(Tetra);


Lamda = [0.25 0.25 0.25];
T = [X(1)-X(4) X(2)-X(4) X(3)-X(4);...
     Y(1)-Y(4) Y(2)-Y(4) Y(3)-Y(4);...
     Z(1)-Z(4) Z(2)-Z(4) Z(3)-Z(4)];
 R = T*Lamda';
 R = R+Vertices(:,4);
 
 r = LMObj('Vertices',R);
 r.Axes = v.RenderAxes;r.Visible = true;