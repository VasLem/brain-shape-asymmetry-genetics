close all; clear all;
addpath('/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/CODE/ANALYSIS/');
addpath('/uz/data/avalok/mic/tmp/SHARED/pclaes4/1KGP/AIDFUNCTIONS/');
studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/';
cd(studypath);
factor = 10000;
neuropath = '/home/pclaes4/Documents/MATLAB/PUBLICGWAS/NEURO/';
%% PARPOOL
try
  parpool('LocalSingle',20);
 catch
end
%% POSSIBLE TABLE HEADERS
listchr = {'CHR' 'chr' 'CHROM'};
listpos = {'BP' 'bp' 'POS' 'pos'};
listP = {'P' 'p' 'PVAL' 'pval' 'Pvalue' 'P-value' 'p-value' 'PVal' 'P_value'};
listSNP = {'SNP' 'snp' 'RSID' 'rsid' 'ID'};
listMarker  = {'variant' 'VARIANT' 'CHR_BP' 'chr_p' 'MarkerName' 'CHR_BP_hg19b37'};
pCrit = 1e-5;
%% NEUROLOGICAL FIRST
neurosubpath = [neuropath 'NEUROLOGICAL/'];
files = dir([neurosubpath '*.tsv']);
nfiles = length(files);
failed = zeros(1,nfiles);
for i=1:nfiles
    try
        %i=1;
        NAME = files(i).name(1:end-4);
        if ~isempty(dir([neurosubpath NAME '.mat'])),continue;end
        disp([ num2str(i) '   ' NAME]);
        Tab = readtable([neurosubpath files(i).name],'FileType','text');
        test = summary(Tab);
        GWAS = [];
        GWAS.TRAIT = NAME;
        GWAS.CHR = [];
        GWAS.POS = [];
        GWAS.P = [];
        % checking for P value info
        list = listP;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               %GWAS.P = uint32(-log10(TMP)*factor);
               GWAS.P = TMP;
               break;
            end
        end
        if iscell(GWAS.P)
           TMP = GWAS.P;
           P = ones(length(TMP),1);
           parfor s=1:length(TMP)
               P(s) = str2double(TMP{s});
           end
           GWAS.P = P;
        end
        % PERFORM AN INITIAL REDUCTION OF THE GWAS DATA
        index = find(GWAS.P<=pCrit);
        GWAS.P = single(-log10(GWAS.P(index)));
        Tab = Tab(index,:);
        % checking for CHR info
        list = listchr;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               GWAS.CHR = TMP;
               break;
            end
        end
        if iscell(GWAS.CHR)
           TMP = GWAS.CHR;
           P = ones(length(TMP),1);
           parfor s=1:length(TMP)
               P(s) = str2double(TMP{s});
           end
           GWAS.CHR = P;
        end
        GWAS.CHR = uint8(GWAS.CHR);
        % checking for POS info
        list = listpos;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               GWAS.POS = TMP;
               break;
            end
        end
        if iscell(GWAS.POS)
           TMP = GWAS.POS;
           P = ones(length(TMP),1);
           parfor s=1:length(TMP)
               P(s) = str2double(TMP{s});
           end
           GWAS.POS = P;
        end
        GWAS.POS = uint32(GWAS.POS);
        % CHECKING IF POS AND CHR ARE NOT EMPTY and otherwise look for
        % variant encoding
        if isempty(GWAS.POS)
            list = listMarker;
            for j=1:1:length(list)
                %j=1;
                if isfield(test,list{j})
                   TMP = eval(['Tab.' list{j}]);
                   nVAR = length(TMP);
                   POS = zeros(nVAR,1,'uint32');
                   CHR = zeros(nVAR,1,'uint8');
                   [path,ID] = setupParForProgress(nVAR);tic;
                   parfor k=1:nVAR
                      str = TMP{k};
                      ind = strfind(str,':');
                      CHR(k) = uint8(str2double(str(1:ind(1)-1)));
                      switch length(ind)
                          case 1
                             POS(k) = uint32(str2double(str(ind(1)+1:end))); 
                          otherwise
                             POS(k) = uint32(str2double(str(ind(1)+1:ind(2)-1)));
                      end
                      parfor_progress;
                   end
                   closeParForProgress(path,ID);toc;
                   GWAS.POS = POS;
                   GWAS.CHR = CHR;
                   break;
                end
            end
         end
         parforSaveToDisk(neurosubpath,[NAME '.mat'],GWAS);
    catch
        failed(i) = 1;
    end
end
%% VOLUMES SECOND
neurosubpath = [neuropath 'VOLUMES/'];
files = dir([neurosubpath '*.tsv']);
nfiles = length(files);
failed = zeros(1,nfiles);
for i=1:nfiles
    try
        %i=8;
        NAME = files(i).name(1:end-4);
        disp([ num2str(i) '   ' NAME]);
        Tab = readtable([neurosubpath files(i).name],'FileType','text');
        test = summary(Tab);
        GWAS = [];
        GWAS.TRAIT = NAME;
        GWAS.CHR = [];
        GWAS.POS = [];
        GWAS.P = [];
        % checking for P value info
        list = listP;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               %GWAS.P = uint32(-log10(TMP)*factor);
               GWAS.P = TMP;
               break;
            end
        end
        
        % PERFORM AN INITIAL REDUCTION OF THE GWAS DATA
        index = find(GWAS.P<=pCrit);
        GWAS.P = single(-log10(GWAS.P(index)));
        Tab = Tab(index,:);
        
        % checking for CHR info
        list = listchr;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               GWAS.CHR = uint8(TMP);
               break;
            end
        end
        % checking for POS info
        list = listpos;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               GWAS.POS = uint32(TMP);
               break;
            end
        end
        % CHECKING IF POS AND CHR ARE NOT EMPTY and otherwise look for
        % variant encoding
        if isempty(GWAS.POS)
            list = listMarker;
            for j=1:1:length(list)
                %j=1;
                if isfield(test,list{j})
                   TMP = eval(['Tab.' list{j}]);
                   nVAR = length(TMP);
                   POS = zeros(nVAR,1,'uint32');
                   CHR = zeros(nVAR,1,'uint8');
                   [path,ID] = setupParForProgress(nVAR);tic;
                   parfor k=1:nVAR
                      str = TMP{k};
                      ind = strfind(str,':');
                      CHR(k) = uint8(str2double(str(1:ind(1)-1)));
                      switch length(ind)
                          case 1
                             POS(k) = uint32(str2double(str(ind(1)+1:end))); 
                          otherwise
                             POS(k) = uint32(str2double(str(ind(1)+1:ind(2)-1)));
                      end
                      parfor_progress;
                   end
                   closeParForProgress(path,ID);toc;
                   GWAS.POS = POS;
                   GWAS.CHR = CHR;
                   break;
                end
            end
         end
         parforSaveToDisk(neurosubpath,[NAME '.mat'],GWAS);
    catch
        failed(i) = 1;
    end
end
%% TXT VOLUMES THIRD
neurosubpath = [neuropath 'VOLUMES/TXT/'];
files = dir([neurosubpath '*.txt']);
nfiles = length(files);
failed = zeros(1,nfiles);
for i=1:nfiles
    try
        %i=2;
        NAME = files(i).name(1:end-4);
        disp([ num2str(i) '   ' NAME]);
        Tab = readtable([neurosubpath files(i).name],'FileType','text');
        test = summary(Tab);
        GWAS = [];
        GWAS.TRAIT = NAME;
        GWAS.CHR = [];
        GWAS.POS = [];
        GWAS.P = [];
        % checking for P value info
        list = listP;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               %GWAS.P = uint32(-log10(TMP)*factor);
               GWAS.P = TMP;
               break;
            end
        end
        
        % PERFORM AN INITIAL REDUCTION OF THE GWAS DATA
        index = find(GWAS.P<=pCrit);
        GWAS.P = single(-log10(GWAS.P(index)));
        Tab = Tab(index,:);
        
        % checking for CHR info
        list = listchr;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               GWAS.CHR = uint8(TMP);
               break;
            end
        end
        % checking for POS info
        list = listpos;
        for j=1:1:length(list)
            %j=1;
            if isfield(test,list{j})
               TMP = eval(['Tab.' list{j}]);
               GWAS.POS = uint32(TMP);
               break;
            end
        end
        % CHECKING IF POS AND CHR ARE NOT EMPTY and otherwise look for
        % variant encoding
        if isempty(GWAS.POS)
            list = listMarker;
            for j=1:1:length(list)
                %j=1;
                if isfield(test,list{j})
                   TMP = eval(['Tab.' list{j}]);
                   nVAR = length(TMP);
                   POS = zeros(nVAR,1,'uint32');
                   CHR = zeros(nVAR,1,'uint8');
                   [path,ID] = setupParForProgress(nVAR);tic;
                   parfor k=1:nVAR
                      str = TMP{k};
                      ind = strfind(str,':');
                      CHR(k) = uint8(str2double(str(1:ind(1)-1)));
                      switch length(ind)
                          case 1
                             POS(k) = uint32(str2double(str(ind(1)+1:end))); 
                          otherwise
                             POS(k) = uint32(str2double(str(ind(1)+1:ind(2)-1)));
                      end
                      parfor_progress;
                   end
                   closeParForProgress(path,ID);toc;
                   GWAS.POS = POS;
                   GWAS.CHR = CHR;
                   break;
                end
            end
         end
         parforSaveToDisk(neurosubpath,[NAME '.mat'],GWAS);
    catch
        failed(i) = 1;
    end
end
%%



%         list = listSNP;
%         for j=1:1:length(list)
%             %j=1;
%             if isfield(test,list{j})
%                TMP = eval(['Tab.' list{j}]);
%                GWAS.RS = TMP;
%                break;
%             end
%         end