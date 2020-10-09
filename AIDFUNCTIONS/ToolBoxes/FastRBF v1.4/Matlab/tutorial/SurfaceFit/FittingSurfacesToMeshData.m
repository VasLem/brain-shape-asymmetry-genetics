% FILENAME: FittingSurfacesToMeshData.m

%cd SurfaceFit;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

Mesh = fastrbf_import('holey_face.obj');
info = fastrbf_checkmesh(Mesh, 'verbose');
f71a = figure_by_name('Figure 7.1(a)');
fastrbf_view(Mesh, 'e');
view(2); colormap([0 0 0]);

f71b = figure_by_name('Figure 7.1(b)');
fastrbf_view(Mesh);
view(2); colormap white; lighting g; camlight left; material shiny;

MeshWithNormals = fastrbf_normalsfrommesh(Mesh);
f72a = figure_by_name('Figure 7.2(a)');
fastrbf_view(MeshWithNormals, 'fv');
view(2); colormap white; lighting g; camlight left; material shiny;

Density = fastrbf_densityfromnormals(MeshWithNormals, 0.5, 5.0);
f72b = figure_by_name('Figure 7.2(b)');
fastrbf_view(Density); view(2);

Density = fastrbf_unique(Density);

rbf = fastrbf_fit(Density, 0.5, 'direct');
f73a = figure_by_name('Figure 7.3(a)');
fastrbf_view(rbf); view(2); colorbar;

rbf = fastrbf_fit(Density, 0.5, 'reduce');
f73b = figure_by_name('Figure 7.3(b)');
fastrbf_view(rbf); view(2); colorbar;

g = fastrbf_grideval(rbf, 0.1, 'spacing', 3);
f74 = figure_by_name('Figure 7.4');
h = slice(g.Value, 40, 27, 15);
set(h, 'edgecolor', 'none', 'facecolor', 'interp')
axis equal; view(3); colorbar;

NewMesh = fastrbf_isosurf(rbf, 2)
f75a = figure_by_name('Figure 7.5(a)');
fastrbf_view(NewMesh, 'e');
view(2); colormap([0 0 0]);

f75b = figure_by_name('Figure 7.5(b)');
fastrbf_view(NewMesh);
view(2); colormap white; lighting g; camlight left; material shiny;

fastrbf_export(NewMesh, 'NewMesh.obj');


% close figures
delete(f71a, f71b, f72a, f72b, f73a, f73b, f74, f75a, f75b);

%cd ..;
