cd('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/DNABIOMETRICS/HAPMAP/HAPMAP/Data');
cd('/Users/dgors0/Documents/MATLAB/Spring 2017/Compare different clusterings')

clear all
close all
clc

load('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/DNABIOMETRICS/HAPMAP/HAPMAP/Data/FACIALDATA/MODULES.mat')
nLevels = 6;

Connectome_meanA_levels = cell(nLevels,1);
Connectome_NmeanA_levels = cell(nLevels,1);
%% Step 1: Calculate the NmeanA and meanA
for level = 1:nLevels,
    fprintf('Calculating level %d ...\n',level);
    nClusters = 2^(level-1);
    Connectome_meanA_levels{level} = zeros(nClusters,nClusters);
    for cl1 = 1:nClusters,
        for cl2 = cl1:nClusters,
            
            % the variable meanA;
            Connectome_meanA_levels{level}(cl1,cl2) = mean(mean(RVMatrix(label(level,:)==cl1,label(level,:)==cl2)));
            Connectome_meanA_levels{level}(cl2,cl1) = Connectome_meanA_levels{level}(cl1,cl2);
        end
    end
    Connectome_NmeanA_levels{level} = Connectome_meanA_levels{level}./sqrt(diag(Connectome_meanA_levels{level})*diag(Connectome_meanA_levels{level})'); 
end

%% Step 2: Plot the matrix version of the Connectome
plotRange = [0 1];
for level = [2:nLevels],
    figure()
    subplot(1,2,1)
    imagesc(Connectome_meanA_levels{level})
    axis square
    caxis(plotRange)
    colorbar
    colormap('jet')
    title({'Connectome: meanA'})
    ylabel(strcat('level',{' '},num2str(level)))
    set(gca,'fontsize',14)
    
    subplot(1,2,2)
    imagesc(Connectome_NmeanA_levels{level})
    axis square
    caxis(plotRange)
    colorbar
    colormap('jet')
    title({'Connectome: NmeanA'})
    set(gca,'fontsize',14)
    set(gcf,'position',[720 850 1250 500]);
end

%% Step 3: Plot the connectome
% find the maximum
MAX = 0;
for level = 2:nLevels,
    temp = triu(Connectome_NmeanA_levels{level},1);
    MAX = max([MAX,diag(Connectome_meanA_levels{level})',temp(:)']);
end

for level = [2:nLevels],
    D_plotConnectome(Connectome_NmeanA_levels{level},diag(Connectome_meanA_levels{level}),[0.5 MAX],[0.5 MAX],level)
end