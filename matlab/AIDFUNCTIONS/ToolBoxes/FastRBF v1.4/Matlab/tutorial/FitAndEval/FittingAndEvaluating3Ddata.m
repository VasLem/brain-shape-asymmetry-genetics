% FILENAME: FittingAndEvaluating3Ddata.m

% Fitting and evaluating 3D data

%cd FitAndEval;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

% test 3D fit, grid eval and export
data = fastrbf_load('data.p3d');
f64 = figure_by_name('Figure 6.4');
colormap(jet(256));
fastrbf_view(data);

data = fastrbf_unique(data);
rbf = fastrbf_fit(data, 0.001);
rbf = fastrbf_fit(data, 0.001, 'reduce', 'verbose');

g = fastrbf_grideval(rbf, 'spacing', 0.02);

[gX, gY, gZ] = fastrbf_gridcoords(g);
gV = permute(g.Value, [2 1 3]);
f65a = figure_by_name('Figure 6.5 (a)');
slice(gY, gX, gZ, gV, .7, .7, .3); colorbar; axis tight;

p = fastrbf_grideval(rbf, 'spacing', 0.02, 'pointlist');
fastrbf_export(p, 'grid.csv', 'txt', 'format', '%x,%y,%z,%v\r\n');

gl = fastrbf_grideval(rbf, 'spacing', 0.02, 'min', [-0.5 -0.5 -0.5], 'max', [1.5 1.5 1.5]);
[glX, glY, glZ] = fastrbf_gridcoords(gl);
glV = permute(gl.Value, [2 1 3]);
f65b = figure_by_name('Figure 6.5 (b)');
slice(glY, glX, glZ, glV, .7, .7, .3); colorbar; axis tight;

% close figures
delete(f64, f65a, f65b);

%cd g..;
