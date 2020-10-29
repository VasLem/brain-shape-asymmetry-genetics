% FILENAME: ImportingAndExportingSurfaceData.m

%cd ImportData;

% Importing and exporting surface data

mesh=fastrbf_import('Hand.obj');

mesh=fastrbf_import('Hand.obj','obj');

fastrbf_export(mesh,'hand.stl');

fastrbf_export(mesh,'hand.stl','stl');

%cd ..;
