%% script to start from folder;

function pcaFROMFolder(Folder,normsize,RefScan,MaskIndex)
% close all;clear all;clear classes;
% Folder = 'C:\MATLAB\WorkTrimmed\Projects\PCAFROMFOLDER\DATA';
% load('C:\MATLAB\WorkTrimmed\Projects\PCAFROMFOLDER\AMv5.mat');
% normsize = true;
iter = 3;
%%
if ~isemtpy(MaskIndex), RefScan = crop(RefScan,'VertexIndex',MaskIndex); end
cd(Folder);
disp('Loading Data');
B = batchCollector;
B.RefScan = clone(RefScan);
B.MaskIndex = MaskIndex;
B.InputFolder = Folder;
%%
process(B);
%%
disp('Computing Model');
Data = B.Data.Shape;
ShapeModel = shapePCA;
ShapeModel.RefScan = clone(RefScan);
getAverage(ShapeModel,Data);
Data = LSGenProcrustes(ShapeModel,Data,normsize,iter);
drawnow;
getModel(ShapeModel,Data);
%% exporting results
mkdir('PCARESULTS');
cd PCARESULTS;
save('Model','ShapeModel');
average = ShapeModel.Average;
save('Average','average');
disp('exporting results');
xlswrite('EigenValues',ShapeModel.EigVal);
xlswrite('EigenVectors',ShapeModel.EigVec);
xlswrite('PCAcoeff',ShapeModel.Tcoeff);
xlswrite('IDs',B.Data.Names');

end