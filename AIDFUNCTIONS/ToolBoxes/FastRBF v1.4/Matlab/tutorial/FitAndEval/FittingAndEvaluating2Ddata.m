% FILENAME: FittingAndEvaluating2Ddata.m

% Fitting and evaluating 2D data

%cd FitAndEval;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

data = fastrbf_import('testdata2d.txt','format','Loc %x %y Tval %v Pval %(*s) \n','comment','%');

data = fastrbf_unique(data);

f61a = figure_by_name('Figure 6.1(a)');
plot(data.Location(1,:), data.Location(2,:), '.');
f61b = figure_by_name('Figure 6.1(b)');
plot3(data.Location(1,:), data.Location(2,:), data.Value, '.');

f61c = figure_by_name('Figure 6.1(c)');
fastrbf_view(data);

rbf = fastrbf_fit(data,0.0001);

rbf = fastrbf_fit(data,0.0001,'verbose');

g = fastrbf_grideval(rbf,'spacing',0.01);

g = fastrbf_grideval(rbf,'spacing',0.01,'accuracy',0.005);

f62a = figure_by_name('Figure 6.2(a)');
fastrbf_view(g);
f62b = figure_by_name('Figure 6.2(b)');
surf(g.Value');

p = fastrbf_grideval(rbf,'spacing',0.01, 'pointlist');
fastrbf_export(p, 'grid.csv', 'txt', 'format', '%x,%y,%v\n');

gl = fastrbf_grideval(rbf, 'spacing', 0.01, 'min', [-0.5 -0.5], 'max', [1.5 1.5])
fastrbf_view(gl);
surf(gl.Value');

g = fastrbf_grideval(rbf,'spacing', [0.01 0.02]);

data = fastrbf_import('testdata2d.txt', 'format', 'Loc %x %y Tval %v Pval %<Pval>v\n', 'comment', '%');
data = fastrbf_unique(data);
rbf = fastrbf_fit(data, 0.0001, 'attr', 'Pval');
g = fastrbf_grideval(rbf, 'spacing', 0.01, 'accuracy', 0.005);

% Iterative fitting
rbf = fastrbf_fit(data, 0.005);
rbf = fastrbf_fit(data, 0.000001, 'initial', rbf);

% Adding new data to an RBF model
% initial fit
rbf = fastrbf_fit(data, 0.0001);
% load new data
newdata = fastrbf_import('newtestdata2d.txt', 'format', 'Loc %x %y Tval %v Pval %<Pval>v\n', 'comment', '%');
data.Location = [data.Location newdata.Location];
data.Value = [data.Value newdata.Value];
data.Pval.Value = [data.Pval.Value newdata.Pval.Value];
% perform new fit
rbf = fastrbf_fit(data, 0.0001, 'initial', rbf);

% close figures
delete(f61a, f61b, f61c, f62a, f62b);

% back to tutorial directory
%cd ..;
