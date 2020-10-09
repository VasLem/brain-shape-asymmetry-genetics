% creating my maps
close all;clear all;
savepath = 'D:\Dropbox\==MATLAB==\myToolboxes\COLORMAPS';
cd(savepath);
%% basic building blocks
cmap = colormap('lines');
linesblock = cmap(1:7,:);
binblock = [0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;0.9 0.9 0.9];

%% single color bars
steps = 64;
linesmaps = nan*zeros(steps,3,size(linesblock,1));
for nc = 1:size(linesblock,1)
    %nc = 1;
    lowcolor = [0 0 0];
    highcolor = linesblock(nc,:);
    map = zeros(steps,3);
    map(1,:) = lowcolor;
    map(end,:) = highcolor;
    for i=1:1:3
        lc = lowcolor(i);
        hc = highcolor(i);
        for s=2:steps-1
            map(s,i) = lc+(s/steps)*(hc-lc);
        end
    end
    linesmaps(:,:,nc) = map;
end

%% opening a test figure

figure;hold on;scatter(rand(1,100),rand(1,100),40,rand(1,100),'filled');
colormap(linesmaps(:,:,1));

binblock = [0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;0.9 0.9 0.9];
steps = 64;
binmaps = nan*zeros(steps,3,size(binblock,1));
for nc = 1:size(binblock,1)
    %nc = 1;
    lowcolor = [0 0 0];
    highcolor = binblock(nc,:);
    map = zeros(steps,3);
    map(1,:) = lowcolor;
    map(end,:) = highcolor;
    for i=1:1:3
        lc = lowcolor(i);
        hc = highcolor(i);
        for s=2:steps-1
            map(s,i) = lc+(s/steps)*(hc-lc);
        end
    end
    binmaps(:,:,nc) = map;
end

figure;hold on;scatter(rand(1,100),rand(1,100),40,rand(1,100),'filled');
colormap(linesmaps(:,:,1));

%% LINES WITH GRAY INSTEAD OF GREEN
mylines = colormap('lines');
index = find(ismember(mylines,[0.4660    0.6740    0.1880],'row'));
mylineswithgray = mylines;
for i=1:1:length(index)
   mylineswithgray(index(i),:) = [0.4 0.4 0.4]; 
end