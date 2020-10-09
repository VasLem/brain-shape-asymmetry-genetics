% FILENAME: Isosurfacing.m

%cd IsoSurf;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

Mesh = fastrbf_import('Hand.obj');
Density = fastrbf_densityfromnormals(Mesh, 0.5, 5.0);
Density = fastrbf_unique(Density);
rbf = fastrbf_fit(Density, 0.25,'reduce');

NewMesh = fastrbf_isosurf(rbf, 1);
f81a = figure_by_name('Figure 8.1(a)');
fastrbf_view(NewMesh, 'fe');
colormap copper; lighting g; camlight left; material shiny; view([20 60]);

NewMesh = fastrbf_isosurf(rbf, 1, 2);

NewMesh = fastrbf_isosurf(rbf, 1, 2, 0.1);

f81b = figure_by_name('Figure 8.1(b)');
NewMesh = fastrbf_isosurf(rbf, 2);
fastrbf_view(NewMesh, 'fe');
colormap copper; lighting g; camlight left; material shiny; view([20 60]);

NewMesh = fastrbf_isosurf(rbf, 1, 'open');
f82a = figure_by_name('Figure 8.2(a)');
fastrbf_view(NewMesh);
colormap copper; lighting g; camlight left; material shiny; view([20 60]);

NewMesh = fastrbf_isosurf(rbf, 1, 'closeplus');
f82b = figure_by_name('Figure 8.2(b)');
fastrbf_view(NewMesh);
colormap copper; lighting g; camlight left; material shiny; view([20 60]);

NewMesh = fastrbf_isosurf(rbf, 1, 'closeminus');
f82c = figure_by_name('Figure 8.2(c)');
fastrbf_view(NewMesh);
colormap copper; lighting g; camlight left; material shiny; view([20 60]);

NewMesh = fastrbf_isosurf(rbf, 2, 'plane');
f83a = figure_by_name('Figure 8.3(a)');
fastrbf_view(NewMesh, 'fe');
colormap copper; lighting g; camlight left; material shiny; view([20 60]);

NewMesh = fastrbf_isosurf(rbf, 2, 'plane');
f83d = figure_by_name('Figure 8.3(d)');
fastrbf_view(NewMesh, 'fe');
view(0, 180); colormap copper; lighting g; camlight left; camzoom(1.5); material shiny; 

rx = [ 1 0 0 0;
0 0 -1 0;
0 1 0 0;
0 0 0 1 ]
NewMesh = fastrbf_isosurf(rbf, 2, 'plane', 'orient', rx);
f83b = figure_by_name('Figure 8.3(b)');
fastrbf_view(NewMesh, 'fe');
colormap copper; lighting g; camlight left; material shiny; view([20 60]);

f83e = figure_by_name('Figure 8.3(e)');
fastrbf_view(NewMesh, 'fe');
view(0, 90); colormap copper; lighting g; camlight left; camzoom(1.5); material shiny;
set(gca, 'XLim', [-23 62], 'YLim', [-20 83]);

NewMesh = fastrbf_isosurf(rbf, 2, 'plane', 'up', [0 1 0]);

rx = [ 1 0 0 0;
0 cos(pi/6) -sin(pi/6) 0;
0 sin(pi/6) cos(pi/6) 0;
0 0 0 1 ];
NewMesh = fastrbf_isosurf(rbf, 2, 'plane', 'orient', rx);
f83c = figure_by_name('Figure 8.3(c)');
fastrbf_view(NewMesh, 'fe');
colormap copper; lighting g; camlight left; material shiny; view(20, 60);

min = [NaN 35 NaN];
NewMesh = fastrbf_isosurf(rbf, 2, 'min', min);
f84a = figure_by_name('Figure 8.4(a)');
fastrbf_view(NewMesh);
colormap copper; lighting g; camlight left; material shiny; view(20, 60);

NewMesh = fastrbf_isosurf(rbf, 2, 'min', min, 'seedlimit', 10);
f84b = figure_by_name('Figure 8.4(b)');
fastrbf_view(NewMesh);
colormap copper; lighting g; camlight left; material shiny; view(20, 60);

NewMesh = fastrbf_isosurf(rbf, 2, 'min', min, 'seedlimit', 0);
f84c = figure_by_name('Figure 8.4(c)');
fastrbf_view(NewMesh);
colormap copper; lighting g; camlight left; material shiny; view(20, 60);

G = fastrbf_grideval(rbf, 'spacing', 2);
mesh = fastrbf_mcubes(G);
f85a = figure_by_name('Figure 8.5(a)');
fastrbf_view(mesh);
colormap copper; lighting g; camlight left; material shiny; view(20, 60);


% close figures
delete(f81a, f81b, f82a, f82b, f82c, f83a, f83b, f83c, f83d, f83e, f84a, f84b, f84c, f85a);

%cd ..;
