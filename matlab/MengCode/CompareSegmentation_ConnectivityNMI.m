clc;clear all;
 
outpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/01_Pheno/Data/';
figpath = '/usr/local/micapollo01/MIC/DATA/STAFF/myuan0/tmp/BrainG2L/01_Pheno/Fig/CompareSegMethods/';

%% Compare Connectivity
ConnectivityTab_DistanceMap = table2array(readtable([outpath 'ConnectivityTab_DistanceMap.xlsx']));
ConnectivityTab_ShapeMap = table2array(readtable([outpath 'ConnectivityTab_ShapeMap.xlsx']));
ConnectivityTab_Smooth = table2array(readtable([outpath 'ConnectivityTab_Smooth.xlsx']));
% coverage 
figure;
x1 = ConnectivityTab_DistanceMap(:,1);
y1 = ConnectivityTab_DistanceMap(:,4);
plot(x1,y1, '-o', 'MarkerSize',10,'LineWidth',3)
hold on
x2 = ConnectivityTab_ShapeMap(:,1);
y2 = ConnectivityTab_ShapeMap(:,4);
plot(x2,y2, '-*', 'MarkerSize',10,'LineWidth',3)
hold on
x3 = ConnectivityTab_Smooth(:,1);
y3= ConnectivityTab_Smooth(:,4);
plot(x3,y3, '-*', 'MarkerSize',10,'LineWidth',3)
legend('%DistanceMatrix','%ShapeMatrix','Smooth Iteration', 'Fontsize', 12,'Location','southeast')
title({'Connectivity:','mean(coverage of the largest cluster within one segment)'},'FontSize', 15);
hold off
saveas(gcf,[figpath 'Connectivity_coverage.jpg']);
% percentage of one-piece segment
figure;
x1 = ConnectivityTab_DistanceMap(:,1);
y1 = ConnectivityTab_DistanceMap(:,3)./ConnectivityTab_DistanceMap(:,2);
plot(x1,y1, '-o', 'MarkerSize',10,'LineWidth',3)
hold on
x2 = ConnectivityTab_ShapeMap(:,1);
y2 = ConnectivityTab_ShapeMap(:,3)./ConnectivityTab_ShapeMap(:,2);
plot(x2,y2, '-*', 'MarkerSize',10,'LineWidth',3)
hold on
x3 = ConnectivityTab_Smooth(:,1);
y3 = ConnectivityTab_Smooth(:,3)./ConnectivityTab_Smooth(:,2);
plot(x3,y3, '-*', 'MarkerSize',10,'LineWidth',3)
legend('%DistanceMatrix','%ShapeMatrix','Smooth Iteration', 'Fontsize', 12,'Location','southeast')
title({'Connectivity:','percentage of one-piece segment'},'FontSize', 15);
hold off
saveas(gcf,[figpath 'Connectivity_onepiece.jpg']);

%% Compare NMI score
NMITab_DistanceMap = table2array(readtable([outpath 'NMITab_DistanceMap.xlsx']));
NMITab_ShapeMap = table2array(readtable([outpath 'NMITab_ShapeMap.xlsx']));
NMITab_Smooth = table2array(readtable([outpath 'NMITab_Smooth.xlsx']));
for i=1:4
figure;
x1 = NMITab_DistanceMap(:,1);
y1 = NMITab_DistanceMap(:,i+2);
plot(x1,y1, '-o', 'MarkerSize',10,'LineWidth',3)
hold on
x2 = NMITab_ShapeMap(:,1);
y2 = NMITab_ShapeMap(:,i+2);
plot(x2,y2, '-*', 'MarkerSize',10,'LineWidth',3)
hold on
x3 = NMITab_Smooth(:,1);
y3 = NMITab_Smooth(:,i+2);
plot(x3,y3, '-x', 'MarkerSize',10,'LineWidth',3)
legend('%DistanceMatrix','%ShapeMatrix','Smooth Iteration', 'Fontsize', 12,'Location','northwest')
title(['NMI score at Level' num2str(i+1)],'FontSize', 15);
hold off
saveas(gcf,[figpath 'NMIscore_Level' num2str(i+1) '.jpg']);
end 
close all;







