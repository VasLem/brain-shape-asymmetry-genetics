% FILENAME: MeshSimplifying.m

%cd Simplify;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

can = fastrbf_load('can.msh');
f87a = figure_by_name('Figure 8.7(a)');
rgb = [can.Red.Value; can.Green.Value; can.Blue.Value];
fastrbf_view(can, 'fe', rgb);

f87b = figure_by_name('Figure 8.7(b)');
fastrbf_view(can, 'f', rgb);

can = fastrbf_normalsfrommesh(can)

density = fastrbf_densityfromnormals(can, 0.2, 5, 'increase', 0.5)
rbf = fastrbf_fit(density, 0.2)
can_rbfmesh = fastrbf_isosurf(rbf, 8)

red = fastrbf_fit(can, 1e-3, 'attr', 'Red');
green = fastrbf_fit(can, 1e-3, 'attr', 'Green');
blue = fastrbf_fit(can, 1e-3, 'attr', 'Blue');

can_rbfmesh = fastrbf_pointeval(red, can_rbfmesh, 'attr', 'Red');
can_rbfmesh = fastrbf_pointeval(green, can_rbfmesh, 'attr', 'Green');
can_rbfmesh = fastrbf_pointeval(blue, can_rbfmesh, 'attr', 'Blue');

rgb = [can_rbfmesh.Red.Value; can_rbfmesh.Green.Value; can_rbfmesh.Blue.Value];

f88a = figure_by_name('Figure 8.8(a)');
fastrbf_view(can_rbfmesh, 'fe', rgb);
f88b = figure_by_name('Figure 8.8(b)');
fastrbf_view(can_rbfmesh, 'f', rgb);

fastrbf_export(can_rbfmesh, rgb, 1, 'can_rbfmesh.wrl');

can_simp_e1 = fastrbf_simplify(can_rbfmesh, rbf, 1);

can_simp_e1 = fastrbf_pointeval(red, can_simp_e1, 'attr', 'Red');
can_simp_e1 = fastrbf_pointeval(green, can_simp_e1, 'attr','Green');
can_simp_e1 = fastrbf_pointeval(blue, can_simp_e1, 'attr', 'Blue');

f89a = figure_by_name('Figure 8.9(a)');
fastrbf_view(can_simp_e1, 'fe', [can_simp_e1.Red.Value; can_simp_e1.Green.Value; can_simp_e1.Blue.Value]);
f89b = figure_by_name('Figure 8.9(b)');
fastrbf_view(can_simp_e1, 'f', [can_simp_e1.Red.Value; can_simp_e1.Green.Value; can_simp_e1.Blue.Value]);

can_simp_f2k = fastrbf_simplify(can_rbfmesh, rbf, -1, 'faces', 2000);

can_simp_f2k = fastrbf_pointeval(red, can_simp_f2k, 'attr', 'Red');
can_simp_f2k = fastrbf_pointeval(green, can_simp_f2k, 'attr','Green');
can_simp_f2k = fastrbf_pointeval(blue, can_simp_f2k, 'attr', 'Blue');

f810a = figure_by_name('Figure 8.10(a)');
fastrbf_view(can_simp_f2k, 'fe', [can_simp_f2k.Red.Value; can_simp_f2k.Green.Value; can_simp_f2k.Blue.Value]);
f810b = figure_by_name('Figure 8.10(b)');
fastrbf_view(can_simp_f2k, 'f', [can_simp_f2k.Red.Value; can_simp_f2k.Green.Value; can_simp_f2k.Blue.Value]);

can_simp_e1_f2k = fastrbf_simplify(can_rbfmesh, rbf, 1, 'faces', 2000);
can_simp_e1_f2k = fastrbf_pointeval(red, can_simp_e1_f2k, 'attr', 'Red');
can_simp_e1_f2k = fastrbf_pointeval(green, can_simp_e1_f2k, 'attr','Green');
can_simp_e1_f2k = fastrbf_pointeval(blue, can_simp_e1_f2k, 'attr', 'Blue');

f_extra_a = figure_by_name('Figure 8.Extra(a)');
fastrbf_view(can_simp_e1_f2k, 'fe', [can_simp_e1_f2k.Red.Value; can_simp_e1_f2k.Green.Value; can_simp_e1_f2k.Blue.Value]);

can_simp = fastrbf_simplifyrgb(can_rbfmesh, rbf, 1, red, green, blue, 0.03);

can_simp = fastrbf_pointeval(red, can_simp, 'attr', 'Red');
can_simp = fastrbf_pointeval(green, can_simp, 'attr','Green');
can_simp = fastrbf_pointeval(blue, can_simp, 'attr', 'Blue');

f811a = figure_by_name('Figure 8.11(a)');
fastrbf_view(can_simp, 'fe', [can_simp.Red.Value; can_simp.Green.Value; can_simp.Blue.Value]);
f811b = figure_by_name('Figure 8.11(b)');
fastrbf_view(can_simp, 'f', [can_simp.Red.Value; can_simp.Green.Value; can_simp.Blue.Value]);

can_simp_p18_c003 = fastrbf_simplifyrgb(can_rbfmesh, rbf, 1e100, red, green, blue, 0.03, 'points', 1800);
f812a = figure_by_name('Figure 8.12(a)');
fastrbf_view(can_simp_p18_c003, 'fe');
f812b = figure_by_name('Figure 8.12(b)');
fastrbf_view(can_simp_p18_c003, 'f');

can_simp = fastrbf_pointeval(red, can_simp, 'attr', 'Red');
can_simp = fastrbf_pointeval(green, can_simp, 'attr', 'Green');
can_simp = fastrbf_pointeval(blue, can_simp, 'attr', 'Blue');

rgb = [can_simp.Red.Value; can_simp.Green.Value; can_simp.Blue.Value];
fastrbf_export(can_simp, rgb, 1, 'can_simp.wrl');

rgb = [can_simp.Red.Value; can_simp.Green.Value; can_simp.Blue.Value];
f_extra_b = figure_by_name('Figure 8.Extra(b)');
fastrbf_view(can_simp, 'f', rgb);

% close figures
delete(f87a, f87b, f88a, f88b, f89a, f89b, f810a, f810b, f_extra_a, f811a, f811b, f812a, f812b, f_extra_b);

%cd ..;
