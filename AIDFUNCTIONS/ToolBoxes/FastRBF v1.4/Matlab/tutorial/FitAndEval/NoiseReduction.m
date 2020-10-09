% FILENAME: NoiseReduction.m

%cd FitAndEval;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

data = fastrbf_load('data.p3d');
data = fastrbf_unique(data);
rbf = fastrbf_fit(data, 0.001, 'reduce');
noise_data = fastrbf_load('noise_.2_data.p3d');
noise_rbf = fastrbf_fit(noise_data, 0.001);
% end preparation

Mesh = fastrbf_isosurf(rbf, 0.03, 0.5);
f68a = figure_by_name('Figure 6.8(a)');
fastrbf_view(Mesh, 'fe');

NoiseMesh = fastrbf_isosurf(noise_rbf, 0.03, 0.5);
f68b = figure_by_name('Figure 6.8(b)');
fastrbf_view(NoiseMesh, 'fe');

SmoothMesh = fastrbf_isosurf(noise_rbf, 0.03, 0.5, 'smooth', 0.1);
f68c = figure_by_name('Figure 6.8(c)');
fastrbf_view(SmoothMesh, 'fe');

rrfit_rbf = fastrbf_fit(noise_data, 0.1, 'errorbar');
RestrictedRangeMesh = fastrbf_isosurf(rrfit_rbf, 0.03, 0.5, 'seedlimit', 0);
f68d = figure_by_name('Figure 6.8(d)');
fastrbf_view(RestrictedRangeMesh, 'fe');

% Low pass filter and 
RestrictedRangeMesh = fastrbf_isosurf(rrfit_rbf, 0.03, 0.5, 'seedlimit', 0, 'smooth', 0.03);

% Asymmetric error bars
data = fastrbf_load('ebardata.p3d');
rbf = fastrbf_fit(data, 0.05, 'errorbar', 'verbose');
evpnts = fastrbf_pointeval(rbf, data);
min(evpnts.Value-data.Value)
max(evpnts.Value-data.Value)

rbf2 = fastrbf_fit(data, 0.02, 'errorbar', 'initial', rbf, 'verbose');
evpnts = fastrbf_pointeval(rbf2, data);
min(evpnts.Value-data.Value)
max(evpnts.Value-data.Value)

rbf3=fastrbf_fit(data, 0.0, 'errorbar', 'initial', rbf2, 'verbose')
evpnts=fastrbf_pointeval(rbf3,data)
min(evpnts.Value-data.Value)
max(evpnts.Value-data.Value)

% Spline smoothing
rhofit_rbf01 = fastrbf_fit(noise_data, 0.1, 'rho', 0.1)
SplineMesh01 = fastrbf_isosurf(rhofit_rbf01, 0.03, 0.5, 'seedlimit', 0);
f69a = figure_by_name('Figure 6.9(a)');
fastrbf_view(SplineMesh01, 'fe');

rhofit_rbf1 = fastrbf_fit(noise_data, 0.1, 'rho', 1.0);
SplineMesh1 = fastrbf_isosurf(rhofit_rbf1, 0.03, 0.5, 'seedlimit', 0);
f69b = figure_by_name('Figure 6.9(b)');
fastrbf_view(SplineMesh1, 'fe');

% Outliers & the confidence parameter
rrfit_rbf = fastrbf_fit(noise_data, 0.1, 'confidence', 0.99, 'errorbar');

% close figures
delete(f68a, f68b, f68c, f68d, f69a, f69b);

%cd ..;
