% FILENAME: ImportingAndExporting3Ddata.m

%cd ImportData;

% Import
data3d=fastrbf_import('testdata3d.txt','format','%x %y %z %(*s) %v %(*s) %a %(*s) %dx %dy %dz %(*s) %<col>dx %<col>dy %<col>dz\n','comment','%')
data3d=fastrbf_import('testdata3d.txt','format','%x %y %z %(*s) %v %(*s) %a %(*s) %dx %dy %dz %(*s) %<R>v %<G>v %<B>v\n','comment','%');

% Export
data = fastrbf_import('testdata3d.txt','format','%x %y %z val %v acc %a grad %dx %dy %dz col %<col>dx %<col>dy %<col>dz\n');
fastrbf_export(data,'newdata3d.txt','txt','format','%x %y %z %(.4f)v %(.4f)a %(.4f)dx %(.4f)dy %(.4f)dz %(.4f)<col>dx %(.4f)<col>dy %(.4f)<col>dz\n');

%cd ..;
