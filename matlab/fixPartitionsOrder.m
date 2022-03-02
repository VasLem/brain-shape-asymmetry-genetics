%fix old clusters to new clusters
% old=load('../results/hierarchicalClusteringDemo/STAGE00DATA/asymmetry_reduction10.bak/levels4/phenotype_varThres80.mat');
% new=load('../results/hierarchicalClusteringDemo/STAGE00DATA/asymmetry_reduction10/levels4/phenotype_varThres80.mat');
% oldSizes = cellfun(@(x)size(x,2), old.clusterPCAPhenoFeatures) + 1;
% newSizes = cellfun(@(x)size(x,2), new.clusterPCAPhenoFeatures);
%%
% oldAvg = cellfun(@(x)mean(x,'all'), old.clusterPCAPhenoFeatures);
% newAvg = cellfun(@(x)mean(x(:,1:end-1),"all"), new.clusterPCAPhenoFeatures);
%%
% A = abs(oldAvg - newAvg');
% M = matchpairs(A,100);
% from = M(:,1);
% to = M(:,2);
% ret=table(from,to);
% writetable(ret,'../part_map.csv', 'Delimiter',' ','WriteRowNames',0,'WriteVariableNames',0);
%%
T = readtable('../part_map.csv',Delimiter=' ');
from = T.Var1;
to = T.Var2;
%%
rootdir = '../results';
% for c=1:length(from)
%     f = from(c);
%     t = to(c);
%     to_rem_str = sprintf('Part%02.f.csv',f);
%     to_add_str = sprintf('Part%02.f.tmp.csv',t);
%     filelist = dir(fullfile(rootdir, ['**/*' sprintf(to_rem_str,f) '*']));
%     for i=1:numel(filelist)
%         oldname = strcat(filelist(i).folder,'/',filelist(i).name);
%         newname = replace(oldname, to_rem_str, to_add_str);
%         movefile(oldname, newname);
%     end
% end
% %%
% for c=1:length(from)
%     f = from(c);
%     t = to(c);
%     to_rem_str = sprintf('par%02.f.log',f);
%     to_add_str = sprintf('par%02.f.tmp.log',t);
%     filelist = dir(fullfile(rootdir, ['**/*' sprintf(to_rem_str,f) '*']));
%     for i=1:numel(filelist)
%         oldname = strcat(filelist(i).folder,'/',filelist(i).name);
%         newname = replace(oldname, to_rem_str, to_add_str);
%         movefile(oldname, newname);
%     end
% end
% %%
% for c=1:length(from)
%     f = from(c);
%     t = to(c);
%     to_rem_str = sprintf('par%02.f.sumstats',f);
%     to_add_str = sprintf('par%02.f.tmp.sumstats',t);
%     filelist = dir(fullfile(rootdir, ['**/*' sprintf(to_rem_str,f) '*']));
%     for i=1:numel(filelist)
%         oldname = strcat(filelist(i).folder,'/',filelist(i).name);
%         newname = replace(oldname, to_rem_str, to_add_str);
%         movefile(oldname, newname);
%     end
% end
%%
% filelist = dir(fullfile(rootdir, ['**/*withPartCCA.mat']));
%%
% for i=11:numel(filelist)
%     disp(i)
%     name = strcat(filelist(i).folder,'/',filelist(i).name);
%     old=load(name);
%     if isfield(old, 'gTLPartStats')
%         oldfield = 'gTLPartStats';
%         oldIfield = 'gTLPartIntStats';
%     else
%         oldfield = 'stats';
%         oldIfield = 'intStats';
%     end
%     if isfield(old.(oldfield), 'chiSqSignificance')
%         oldsigfield = 'chiSqSignificance';
%     else
%         oldsigfield = 'chisqSignificance';
%     end
%     new = cell(2,1);
%     for k=1:2
%         if k ==1
%             oldstr = old.(oldfield);
%         else
%             oldstr = old.(oldIfield);
%         end
%         oldstrsig = oldstr.(oldsigfield);
%         new{k}.chisqSignificance = zeros(size(oldstrsig));
%         new{k}.chisq = zeros(size(oldstr.chisq));
%         new{k}.df = zeros(size(oldstr.df));
%         for c=1:length(from)
%             f = from(c);
%             t = to(c);
%             oldstrsig =  oldstr.(oldsigfield);
%             new{k}.chisqSignificance(:,t) = oldstrsig(:,f);
%             new{k}.chisq(:,t) = oldstr.chisq(:,f);
%             new{k}.df(:,t) = oldstr.df(:,f);
%         end
%     end
%     stats = new{1};
%     intStats = new{2};
%     save(name,'stats','intStats');
% end
