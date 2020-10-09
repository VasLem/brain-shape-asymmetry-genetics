% Testing Toolbox Graph

cd C:\MATLAB\Work\toolbox_graph\off;

plot_mesh(vertex, faces);
name = 'elephant-50kv';
options.name = name; % useful for displaying
[vertex,faces] = read_mesh(name);

Oli = meshObj('Vertices',vertex,'Faces',faces);

%% 
load('C:\MATLAB\Work\RCHNormalAverage\DATA\10mF_70.mat');
load('C:\MATLAB\Work\RCHNormalAverage\DATA\MAPPED\10mF_56.mat')
vertex = Scan.Vertices;
faces = Scan.Faces;

options.curvature_smoothing = 10;
options.verb = 0;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,faces,options);

viewer(Scan)
Scan.Value = Cmin';
Scan.ColorMode = 'Indexed';

Scan = clone(Scan);
viewer(Scan)
Scan.Value = Cmax;
Scan.ColorMode = 'Indexed';

Scan = clone(Scan);
viewer(Scan)
Scan.Value = Cmean;
Scan.ColorMode = 'Indexed';

Scan = clone(Scan);
viewer(Scan)
Scan.Value = perform_saturation(Cgauss,1.2);
Scan.ColorMode = 'Indexed';


Scan = clone(Scan);
viewer(Scan)
Scan.Value = perform_saturation(abs(Cmin)+abs(Cmax),1.2);
Scan.ColorMode = 'Indexed';

%% Detecting Nose Tip

Tk = 0.0025;
Th = 0.00005;

indexK = find(Cgauss<Tk);
indexH = find(Cmean>Th);
indexNoseTip = intersect(indexK,indexH);

NoseTip = LMObj('Vertices',Scan.Vertices(:,indexNoseTip));

%% Parameterization

load('C:\MATLAB\Work\RCHNormalAverage\DATA\MAPPED\10mF_56.mat')
vertex = Scan.Vertices;
faces = Scan.Faces;

A = triangulation2adjacency(faces);
% you can try with other boundary type
options.boundary = 'circle';
% you can try with other Laplacians
options.laplacian = 'distance';
%The fixed boundary parameterization is guaranteed to be valid is the mesh is homeomorphic to a disk and if the boundary of the planar domain is convex. This is the Tutte embedding theorem. 

% compute the layout in 2D
options.method = 'parameterization';
options.verb = 0;
xy = compute_parameterization(vertex,faces,options);
% display the parameterization
clf;
subplot(1,2,1);
plot_mesh(vertex,faces,options); shading faceted; axis tight;
subplot(1,2,2);
plot_graph(A,xy,'k.-'); axis tight;
title('Parameterization');


