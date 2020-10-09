% FILENAME: Anti_Aliasing.m

% Fitting and evaluating noisy data

%cd FitAndEval;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

data = fastrbf_load('data.p3d');
data = fastrbf_unique(data);
rbf = fastrbf_fit(data, 0.001, 'reduce');
% end preparation

noise_data = fastrbf_load('noise_.2_data.p3d');

f66a = figure_by_name('Figure 6.6 (a)');
colormap(jet(256));
fastrbf_view(noise_data);

noise_rbf = fastrbf_fit(noise_data, 0.001);

g = fastrbf_grideval(noise_rbf, 'spacing', 0.02)
% display
[gX, gY, gZ] = fastrbf_gridcoords(g);
gV = permute(g.Value, [2 1 3]);
f66b = figure_by_name('Figure 6.6 (b)');
slice(gY, gX, gZ, gV, .7, .7, .3); colorbar; axis tight;

g = fastrbf_grideval(noise_rbf, 'spacing', 0.1)
% display
[gX, gY, gZ] = fastrbf_gridcoords(g);
gV = permute(g.Value, [2 1 3]);
f67a = figure_by_name('Figure 6.7 (a)');
slice(gY, gX, gZ, gV, .7, .7, .3); colorbar; axis tight;

g = fastrbf_grideval(noise_rbf, 'spacing', 0.1, 'smooth', 0.1)
% display
[gX, gY, gZ] = fastrbf_gridcoords(g);
gV = permute(g.Value, [2 1 3]);
f67b = figure_by_name('Figure 6.7 (b)');
slice(gY, gX, gZ, gV, .7, .7, .3); colorbar; axis tight;

g = fastrbf_grideval(rbf, 'spacing', 0.1)
% display
[gX, gY, gZ] = fastrbf_gridcoords(g);
gV = permute(g.Value, [2 1 3]);
f67c = figure_by_name('Figure 6.7 (c)');
slice(gY, gX, gZ, gV, 1, 1, .3); colorbar; axis tight;

% close figures
delete(f66a, f66b, f67a, f67b, f67c);

%cd ..;
