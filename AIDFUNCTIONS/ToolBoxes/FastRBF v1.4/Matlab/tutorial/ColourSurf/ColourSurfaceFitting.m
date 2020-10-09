% FILENAME: ColourSurfaceFitting.m

%cd ColourSurf;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

weta_flat = imread('weta.png');
f76a = figure_by_name('Figure 7.6(a)');
image(weta_flat);

weta_dome = fastrbf_import('weta_dome.obj');

% Hint text
rgb = [weta_dome.Red.Value; weta_dome.Green.Value; weta_dome.Blue.Value]/255;
f76b = figure_by_name('Figure 7.6(b)');
fastrbf_view(weta_dome, 'f', rgb);
view(90,0);

weta_dome = fastrbf_normalsfrommesh(weta_dome);
density=fastrbf_densityfromnormals(weta_dome, 0.001, 0.05);
weta.rbf = fastrbf_fit(density, 1e-3,'reduce');
weta.mesh = fastrbf_isosurf(weta.rbf, 0.01);

weta.Red =fastrbf_fit(weta_dome,0.3,'direct','attr', 'Red');
weta.Green=fastrbf_fit(weta_dome,0.3,'direct','attr', 'Green');
weta.Blue=fastrbf_fit(weta_dome,0.3,'direct','attr', 'Blue');

weta.mesh = fastrbf_pointeval(weta.Red, weta.mesh, 0.1, 'attr', 'Red');
weta.mesh = fastrbf_pointeval(weta.Green, weta.mesh, 0.1, 'attr', 'Green');
weta.mesh = fastrbf_pointeval(weta.Blue, weta.mesh, 0.1, 'attr', 'Blue');

rgb = [weta.mesh.Red.Value; weta.mesh.Green.Value; weta.mesh.Blue.Value]/255;

f77a = figure_by_name('Figure 7.7(a)');
fastrbf_view(weta.mesh, 'f', rgb);
view(90, 0);

weta.iso = fastrbf_isosurf(weta.rbf, 0.01);
weta.iso = fastrbf_pointeval(weta.Red, weta.iso, 0.1, 'attr', 'Red');
weta.iso = fastrbf_pointeval(weta.Green, weta.iso, 0.1, 'attr', 'Green');
weta.iso = fastrbf_pointeval(weta.Blue, weta.iso, 0.1, 'attr', 'Blue');

rgb = [weta.iso.Red.Value; weta.iso.Green.Value; weta.iso.Blue.Value]/255;
f77b = figure_by_name('Figure 7.7(b)');
fastrbf_view(weta.iso, 'f', rgb);
view([90 0]);

% close figures
delete(f76a, f76b, f77a, f77b);

%cd ..;

