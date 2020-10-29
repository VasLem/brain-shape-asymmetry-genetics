% FILENAME: ImportingAndExportingFormatedTextData.m

%cd ImportData;

% preparation
% give access to figure_by_name function
tmp_p = pwd; cd ..; addpath(pwd); cd(tmp_p);

% Read formated text data

data2d = fastrbf_import('testdata2d.txt', 'format', '%(*s) %x %y %(*s) %<Tvalue>v %(*s) %<Pvalue>v\n', 'comment', '%');

data2d = fastrbf_import('testdata2d.txt', 'format', 'Loc %x %y %(*s) %<Tvalue>v %(*s) %<Pvalue>v\n');

data2d = fastrbf_import('testdata2d.txt', 'format', 'Loc %x %y Tval %<Tvalue>v Pval %<Pvalue>v\n');

data2d = fastrbf_import('testdata2d.txt','format','Loc %x %y Tval %<Tvalue>v');

data2d = fastrbf_import('testdata2d.txt','format','Loc %x %y Tval %(*s) Pval %<Pvalue>v\n');

%data=fastrbf_import('testdata2d.txt','format','Loc %x %y Tval %<Tvalue>v Pval %<Pvalue>v\n')
data = fastrbf_import('testdata2d.txt','format','Loc %x %y Tval %v Pval %<Pvalue>v\n')
f51a = figure_by_name('Figure 5.1(a)');
plot3(data.Location(1,:), data.Location(2,:), data.Value, '.');

f51b = figure_by_name('Figure 5.1(b)');
plot3(data.Location(1,:), data.Location(2,:), data.Pvalue.Value, '.');

data = fastrbf_unique(data)
%Save data in FastRBF native format
fastrbf_save(data, 'tut2ddata')
%Fit RBF to Tval & evaluate on grid
rbf_Tval = fastrbf_fit(data,0.0001)
Tgrid = fastrbf_grideval(rbf_Tval, 'spacing', 0.01)
%Display interpolated data
[x y] = fastrbf_gridcoords(Tgrid);
f51c = figure_by_name('Figure 5.1(c)');
surf(x, y, Tgrid.Value, 'facecolor', 'interp');
%Fit RBF to Pval & evaluate on grid
rbf_Pval = fastrbf_fit(data,0.0001,'attr','Pvalue')
Pgrid = fastrbf_grideval(rbf_Pval,'spacing',0.01)
%Display interpolated data
[x y] = fastrbf_gridcoords(Pgrid);
f51d = figure_by_name('Figure 5.1(d)');
surf(x, y, Pgrid.Value, 'facecolor', 'interp');

[x y] = fastrbf_gridvec(Tgrid);
[xx, yy]=ndgrid(x,y);
Pgrid_plist.Location=[xx(:), yy(:)]';
Tgrid_plist.Location=Pgrid_plist.Location;
Tgrid_plist.Value=Tgrid.Value(:)';
Pgrid_plist.Value=Pgrid.Value(:)';
%Export point list data as formatted text
fastrbf_export(Tgrid_plist,'Tgrid.txt','txt','format','%x %y %v\n')
fastrbf_export(Pgrid_plist,'Pgrid.txt','txt','format','%x %y %v\n')

[x y] = fastrbf_gridvec(Tgrid);
[xx, yy]=ndgrid(x,y);
TPgrid_plist.Location=[xx(:), yy(:)]';
TPgrid_plist.Tvalue.Value=Tgrid.Value(:)';
TPgrid_plist.Pvalue.Value=Pgrid.Value(:)';
%Export point list data as formatted text
fastrbf_export(TPgrid_plist,'TPgrid.txt','txt','format','%x %y %<Tvalue>v %<Pvalue>v\n');

%Generate grid data as a point list
Tgrid_plist=fastrbf_grideval(rbf_Tval,'pointlist','spacing', 0.01);
Pgrid_plist=fastrbf_grideval(rbf_Pval,'pointlist','spacing', 0.01);
%Export point list data as formatted text
fastrbf_export(Tgrid_plist,'Tgrid.txt','txt','format','%x %y %v\n');
fastrbf_export(Pgrid_plist,'Pgrid.txt','txt','format','%x %y %v\n');

data=fastrbf_import('testdata2d.txt','format','Loc %x %y Tval %<Tvalue>v Pval %<Pvalue>v\n');
fastrbf_export(data,'newdata.txt','txt','format','%x %y %(.4f)<Tvalue>v %(.4f)<Pvalue>v\n');

% close figures
delete(f51a, f51b, f51c, f51d);

%cd ..;
