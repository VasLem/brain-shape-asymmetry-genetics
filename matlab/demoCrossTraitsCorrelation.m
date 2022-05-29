%Aggregate all munged sumstats into one big mat file based on the available
%SNPs in 1000G

labels ={'OTHER_TRAITS_GWAS','OTHER_ASYMMETRY_GWAS','BRAIN_SHAPE_CHR'};
for lInd=1:1%:length(labels)
    lab = labels{lInd};
    INPUT_DIR = ['../SAMPLE_DATA/' lab];
    disp(['Label:',lab])
    if lInd ~=3 
        
        allfiles = dir(INPUT_DIR);
        % Get a logical vector that tells which is a directory.
        dirFlags = [allfiles.isdir];
        % Extract only those that are directories.
        subFolders = allfiles(dirFlags); % A structure with extra info.
        % Get only the folder names into a cell array.
        subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..
    else
        allfiles = dir(INPUT_DIR);
        subFolderNames = {allfiles(3:end).name};
    end

    SYN = struct;
    cnt = 1;
    for traitInd=1:length(subFolderNames)
        trait = subFolderNames{traitInd};
        
        if lInd ~= 3
            inpstatsFile = [INPUT_DIR '/' trait '/munged.sumstats.gz'];
        else
            inpstatsFile = [INPUT_DIR '/' trait];
            trait = trait(1:end-7);
        end
        
        disp(['Trait:', trait])
        if ~isfile(inpstatsFile), disp([inpstatsFile, ' does not exist']); continue; end
        unzippedFile = gunzip(inpstatsFile);
        if lInd ~= 3
            x = readtable(unzippedFile{1},FileType="text");
        else
            x = readtable(unzippedFile{1},FileType="text",Delimiter=',');
        end
        delete(unzippedFile{1});

        if lInd ~= 3
            p = x.P;
            if traitInd == 1
                SYN.RS = x.SNP;
                SYN.P = p;
            else
                [~,ind21] = vlookupFast(x.SNP, SYN.RS);
                SYN.P = cat( 2, SYN.P, p(ind21));
            end
        else
            if traitInd == 1
                SYN.RS = x.rsID;
                SYN.P = x.P_value;
            else
                [~,ind21] = vlookupFast(x.rsID, SYN.RS);
                SYN.P = cat(2, SYN.P, x.P_value(ind21));
            end
        end
        SYN.NAMES{cnt} = trait;
        cnt=cnt+1;
    end
    out_dir = '../results/correlation_sources/';
    if ~isfolder(out_dir), mkdir(out_dir); end
    save([out_dir lab '.mat'], 'SYN', '-v7.3');
end
%%
%%

