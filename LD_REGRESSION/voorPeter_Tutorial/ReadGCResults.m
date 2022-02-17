studypath = '/usr/local/avalok/tmp/hhoske1/sumstats/CCA/';
studypath2 = '/usr/local/avalok/tmp/hhoske1/sumstats/tmprg/pub/';
cd(studypath)

path1 = [studypath 'US/results/face/'];
path2 = [studypath 'META/results/face/'];
path3 = [studypath 'US/results/public/'];
path4 = [studypath 'META/results/public/'];
nMOD = 63;
%% face
h2_T1 = nan(nMOD,nMOD);
h2_SE_T1 = nan(nMOD,nMOD);
h2_T2 = nan(nMOD,nMOD);
h2_SE_T2 = nan(nMOD,nMOD);
GC = nan(nMOD,nMOD);
GC_SE = nan(nMOD,nMOD);
GC_Z = nan(nMOD,nMOD);
GC_P = nan(nMOD,nMOD);

path = path1;
% path = path2;

remind = [];
for i=1:nMOD-1
    disp(i)
    for j=i+1:nMOD
        in = readcell([path 'GC_' num2str(i) '_' num2str(j) '.log'],'FileType','text');
        
        idx = find(contains(in(:,1),'Heritability of phenotype 1'));
        tmp = in(idx+2,:);
        tmp = tmp(cellfun(@ischar,tmp));% Remove missings
        tmp = cellfun(@(x) split(x,' '),tmp,'UniformOutput',false);
        tmp = vertcat(tmp{:});                    
        h2_T1(i,j) = str2double(tmp(end-1));
        h2_SE_T1(i,j) = str2double(erase(tmp(end),{'(' ')'}));
        
        idx = find(contains(in(:,1),'Heritability of phenotype 2/2'));
        tmp = in(idx+2,:);
        tmp = tmp(cellfun(@ischar,tmp));% Remove missings
        tmp = cellfun(@(x) split(x,' '),tmp,'UniformOutput',false);
        tmp = vertcat(tmp{:});                    
        h2_T2(i,j) = str2double(tmp(end-1));
        h2_SE_T2(i,j) = str2double(erase(tmp(end),{'(' ')'}));
                   
        idx = find(contains(in(:,1),'Summary of Genetic Correlation Results'));  
        tmp = in(idx+1:idx+2,:);
        tmp = tmp(cellfun(@ischar,tmp));% Remove missings
        tmp = split(tmp);
        idx2 = find(contains(tmp(1,:),'rg'));
        GC(i,j) = str2double(tmp(2,idx2));
        GC_SE(i,j) = str2double(tmp(2,idx2+1));
        GC_Z(i,j) = str2double(tmp(2,idx2+2));
        GC_P(i,j) = str2double(tmp(2,idx2+3));   
        
        tmp = in(:);
        tmp = tmp(cellfun(@ischar,tmp));
        if sum(contains(tmp,'rg out of bounds')), remind = [remind; i j]; end
        
    end
end
% isequal(h2_T1(2:end-1,end),h2_T2(1,2:end-1)')
h2 = h2_T1(:,end);
h2(end) = h2_T2(1,end);
% isequal(h2_SE_T1(2:end-1,end),h2_SE_T2(1,2:end-1)')
h2_SE = h2_SE_T1(:,end);
h2_SE(end) = h2_SE_T2(1,end);

switch path
    case path1
        save([studypath 'US/results/GC_US_FACES.mat'],'h2','h2_SE','GC','GC_SE','GC_Z','GC_P','remind')
    case path2
        save([studypath 'META/results/GC_META_FACES.mat'],'h2','h2_SE','GC','GC_SE','GC_Z','GC_P','remind')
end
%% public
files = dir([studypath2 '*.txt']);
filenames = split([files(:).name],'.txt');
filenames = filenames(2:end-1);
nPE = length(filenames);

in = readtable([studypath 'GC_PUBGWAS_INFO.xlsx']);
tmpin = table2array(in(:,1:2));
tmpin(:,1) = cellfun(@lower,tmpin(:,1),'UniformOutput',false);
tmppheno = join(tmpin,'_');
[m,indm] = ismember(filenames,tmppheno);
pheno = in(indm,:);

h2_face = nan(nMOD,nPE);
h2_SE_face = nan(nMOD,nPE);
h2_public = nan(nMOD,nPE);
h2_SE_public = nan(nMOD,nPE);
GC = nan(nMOD,nPE);
GC_SE = nan(nMOD,nPE);
GC_Z = nan(nMOD,nPE);
GC_P = nan(nMOD,nPE);

% path = path3;
path = path4;

remind = [];
for i=1:nMOD
    disp(i)
    for j=1:nPE
        in = readcell([path 'GC_' num2str(i) '_' filenames{j} '.txt.log'],'FileType','text');

        if sum(contains(in(:,1),'ERROR computing rg')), remind = [remind; num2cell(i) filenames(j)]; continue; end
        
        idx = find(contains(in(:,1),'Heritability of phenotype 1'));
        tmp = in(idx+2,:);
        tmp = tmp(cellfun(@ischar,tmp));% Remove missings
        tmp = cellfun(@(x) split(x,' '),tmp,'UniformOutput',false);
        tmp = vertcat(tmp{:});                    
        h2_face(i,j) = str2double(tmp(end-1));
        h2_SE_face(i,j) = str2double(erase(tmp(end),{'(' ')'}));
        
        idx = find(contains(in(:,1),'Heritability of phenotype 2/2'));
        tmp = in(idx+2,:);
        tmp = tmp(cellfun(@ischar,tmp));% Remove missings
        tmp = cellfun(@(x) split(x,' '),tmp,'UniformOutput',false);
        tmp = vertcat(tmp{:});                    
        h2_public(i,j) = str2double(tmp(end-1));
        h2_SE_public(i,j) = str2double(erase(tmp(end),{'(' ')'}));
        
        idx = find(contains(in(:,1),'Summary of Genetic Correlation Results'));  
        tmp = in(idx+1:idx+2,:);
        tmp = tmp(cellfun(@ischar,tmp));% Remove missings
        tmp = split(tmp);
        idx2 = find(contains(tmp(1,:),'rg'));
        GC(i,j) = str2double(tmp(2,idx2));
        GC_SE(i,j) = str2double(tmp(2,idx2+1));
        GC_Z(i,j) = str2double(tmp(2,idx2+2));
        GC_P(i,j) = str2double(tmp(2,idx2+3));   
        
        tmp = in(:);
        tmp = tmp(cellfun(@ischar,tmp));
        if sum(contains(tmp,'rg out of bounds')), remind = [remind; num2cell(i) filenames(j)]; end                
    end
end
h2_face = mean(h2_face,2);
h2_SE_face = mean(h2_SE_face,2);
h2_public = mean(h2_public);
h2_SE_public = mean(h2_SE_public);

switch path
    case path3
         save([studypath 'US/results/GC_US_PUBLIC.mat'],'h2_face','h2_SE_face','h2_public','h2_SE_public','GC','GC_SE','GC_Z','GC_P','remind','pheno')
    case path4
        save([studypath 'META/results/GC_META_PUBLIC.mat'],'h2_face','h2_SE_face','h2_public','h2_SE_public','GC','GC_SE','GC_Z','GC_P','remind','pheno')
end
%%
