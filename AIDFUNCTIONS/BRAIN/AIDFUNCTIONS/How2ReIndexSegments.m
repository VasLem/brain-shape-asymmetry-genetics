% learn to reindex segments according to quadrant
close all;clear all;
studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/'
cd(studypath);
%% FACIAL SEGMENTS
close all;clear all;
studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/'
cd(studypath);
HI = HierarchicalInterface;
ReOrg.Index = [1 2 3 4 sort(getAllChildren(HI,4)) 5 sort(getAllChildren(HI,5)) 6 sort(getAllChildren(HI,6)) 7 sort(getAllChildren(HI,7))];
ReOrg.Label = [1 2 3 4*ones(1,length([4 sort(getAllChildren(HI,4))])) 5*ones(1,length([5 sort(getAllChildren(HI,5))])) 6*ones(1,length([6 sort(getAllChildren(HI,6))])) 7*ones(1,length([7 sort(getAllChildren(HI,7))]))];
FaceReOrg = ReOrg;
save('FaceReOrgInfo.mat','FaceReOrg');
%% LEFT HEMISPHERE
close all;clear all;
studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/'
cd(studypath);
Segments = load('./LH/RENDERMATERIAL.mat');
ReOrg.Index = [1 2 3 4 sort(getAllChildren(Segments.UHI,4)) 5 sort(getAllChildren(Segments.UHI,5)) 6 sort(getAllChildren(Segments.UHI,6)) 7 sort(getAllChildren(Segments.UHI,7))];
ReOrg.Label = [1 2 3 4*ones(1,length([4 sort(getAllChildren(Segments.UHI,4))])) 5*ones(1,length([5 sort(getAllChildren(Segments.UHI,5))])) 6*ones(1,length([6 sort(getAllChildren(Segments.UHI,6))])) 7*ones(1,length([7 sort(getAllChildren(Segments.UHI,7))]))];
index = find(Segments.UMASK);
ReOrg.SegIndex = 0*ReOrg.Index;
for i=1:length(index)
    % i=1;
    ReOrg.SegIndex(ReOrg.Index==index(i)) = i;
end
ind = find(ReOrg.SegIndex);
ReOrg.SegIndex = ReOrg.SegIndex(ind);
ReOrg.SegLabel = ReOrg.Label(ind);
save('./LH/ReOrgInfo.mat','ReOrg');
%% LEFT HEMISPHERE
close all;clear all;
studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/'
cd(studypath);
Segments = load('./RH/RENDERMATERIAL.mat');
ReOrg.Index = [1 2 3 4 sort(getAllChildren(Segments.UHI,4)) 5 sort(getAllChildren(Segments.UHI,5)) 6 sort(getAllChildren(Segments.UHI,6)) 7 sort(getAllChildren(Segments.UHI,7))];
ReOrg.Label = [1 2 3 4*ones(1,length([4 sort(getAllChildren(Segments.UHI,4))])) 5*ones(1,length([5 sort(getAllChildren(Segments.UHI,5))])) 6*ones(1,length([6 sort(getAllChildren(Segments.UHI,6))])) 7*ones(1,length([7 sort(getAllChildren(Segments.UHI,7))]))];
index = find(Segments.UMASK);
ReOrg.SegIndex = 0*ReOrg.Index;
for i=1:length(index)
    % i=1;
    ReOrg.SegIndex(ReOrg.Index==index(i)) = i;
end
ind = find(ReOrg.SegIndex);
ReOrg.SegIndex = ReOrg.SegIndex(ind);
ReOrg.SegLabel = ReOrg.Label(ind);
save('./RH/ReOrgInfo.mat','ReOrg');
%%
close all;clear all;
studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/PHENOTYPING/'
cd(studypath);
Segments = load('./SWH/RENDERMATERIAL.mat');
ReOrg.Index = [1 2 3 4 sort(getAllChildren(Segments.UHI,4)) 5 sort(getAllChildren(Segments.UHI,5)) 6 sort(getAllChildren(Segments.UHI,6)) 7 sort(getAllChildren(Segments.UHI,7))];
ReOrg.Label = [1 2 3 4*ones(1,length([4 sort(getAllChildren(Segments.UHI,4))])) 5*ones(1,length([5 sort(getAllChildren(Segments.UHI,5))])) 6*ones(1,length([6 sort(getAllChildren(Segments.UHI,6))])) 7*ones(1,length([7 sort(getAllChildren(Segments.UHI,7))]))];
index = find(Segments.UMASK);
ReOrg.SegIndex = 0*ReOrg.Index;
ReOrg.SegLabel = 0*ReOrg.Index;
for i=1:length(index)
    % i=7;
    ind = find(ReOrg.Index==index(i));
    ReOrg.SegIndex(ind) = i;
    ReOrg.SegLabel(ind) = ReOrg.Label(ind);
end
ind = find(ReOrg.SegIndex);
ReOrg.OrderedSegIndex = ReOrg.SegIndex(ind);
ReOrg.OrderedSegLabel = ReOrg.Label(ind);
nSegments = length(ReOrg.OrderedSegIndex);
ReOrg.OriginalSegIndex = 1:nSegments;
ReOrg.OriginalSegLabel = zeros(1,nSegments);
for i=1:nSegments
    ind = find(ReOrg.OrderedSegIndex == ReOrg.OriginalSegIndex(i));
    ReOrg.OriginalSegLabel(i) = ReOrg.OrderedSegLabel(ind);
end
save('./SWH/ReOrgInfo.mat','ReOrg');
%%
index = find(Segments.UMASK);
ReOrg.SegIndex = 0*ReOrg.Index;
for i=1:length(index)
    % i=1;
    ReOrg.SegIndex(ReOrg.Index==index(i)) = i;
end
ind = find(ReOrg.SegIndex);
ReOrg.SegIndex = ReOrg.SegIndex(ind);
ReOrg.SegLabel = ReOrg.Label(index);




