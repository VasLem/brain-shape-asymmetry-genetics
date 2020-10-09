% FILENAME: SurfacingNoisyData.m

%cd Smoothing;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

Scan = fastrbf_import('Fountain.ptx')
f78 = figure_by_name('Figure 7.8');
fastrbf_view(Scan);
view(2);

Scan = fastrbf_crop(Scan, [0.100 -0.900 NaN], [0.800 -0.200 NaN]);
f79 = figure_by_name('Figure 7.9');
fastrbf_view(Scan);
view(2); colorbar;

Scan = fastrbf_normalsfromscan(Scan);

Density = fastrbf_densityfromnormals(Scan, 0.006,0.030);

Density = fastrbf_unique(Density);

Rbf1 = fastrbf_fit(Density, 0.002);
Mesh = fastrbf_isosurf(Rbf1, 0.005);
f710a = figure_by_name('Figure 7.10(a)');
fastrbf_view(Mesh); view(2);
colormap gray; lighting g; camlight right; material shiny;

Mesh1 = fastrbf_isosurf(Rbf1, 0.005, 'smooth', 0.02);
f710b = figure_by_name('Figure 7.10(b)');
fastrbf_view(Mesh1); view(2);
colormap gray; lighting g; camlight right; material shiny;

Rbf2 = fastrbf_fit(Density, 0.004, 'errorbar');
Mesh2 = fastrbf_isosurf(Rbf2, 0.005);
f711a = figure_by_name('Figure 7.11(a)');
fastrbf_view(Mesh2); view(2);
colormap gray; lighting g; camlight right; material shiny;

Rbf3 = fastrbf_fit(Density, 0.002, 'rho', 0.1);
Mesh3 = fastrbf_isosurf(Rbf3, 0.005);
f711b = figure_by_name('Figure 7.11(b)');
fastrbf_view(Mesh3); view(2);
colormap gray; lighting g; camlight right; material shiny;

Rbf3 = fastrbf_fit(Density, 0.003, 'errorbar', 'initial', Rbf2);
Mesh3 = fastrbf_isosurf(Rbf3, 0.005);
f711c = figure_by_name('Figure 7.11(c)');
fastrbf_view(Mesh3); view(2);
colormap gray; lighting g; camlight right; material shiny;

% close figures
delete(f78, f79, f710a, f710b, f711a, f711b, f711c);

%cd ..;
