%% Setting Paths
close all;clear all;clear classes;
InputPath = 'D:\==PCLAES==\Work\==Papers==\MyPapers\==3DSurgeryEvaluationPlanning==\IEEEDysmorphometry\Biostatistics\FINALIMAGESBMP';
OutputPath = 'D:\==PCLAES==\Work\==Papers==\MyPapers\==3DSurgeryEvaluationPlanning==\IEEEDysmorphometry\Biostatistics\FINALIMAGESEPS';
cd(InputPath)
Files = dir('*.*');
nrFiles = length(Files);
%%
for i = 3:1:nrFiles
    cd(InputPath);
    %im = imread(Files(i).name);
    imshow(Files(i).name);
    cd(OutputPath);
    print(gcf,'-depsc','-r300',Files(i).name(1:end-4));
    close all;
end
