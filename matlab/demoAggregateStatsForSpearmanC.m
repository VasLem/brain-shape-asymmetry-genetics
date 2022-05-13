clear all
addpath(genpath('AIDFUNCTIONS'));

%%
datasets = {'STAGE00DATA', 'BATCH2_2021_DATA','joinedDatasets'};
for datasetInd=1:3
    dataset = datasets{datasetInd};
    INPUT_DIR = ['../results/asymmetry/meta_analysis/' dataset '/mean_imputed/not_subsampled/'];
    PAR = struct;
    for part=1:31
        disp(['Partition:', num2str(part)]);
        inpStatsFile = [INPUT_DIR 'CCAPart' sprintf('%02d', part), '.csv.gz'];
        unzippedFile = gunzip(inpStatsFile);
        x = readtable(unzippedFile{1},FileType="text",Delimiter=',');
        delete(unzippedFile{1});
        if part == 1
            PAR.RS = x.rsID;
            try
                PAR.POS = x.position;
            catch
            end
            PAR.CHR = x.chromosome;
            PAR.P = x.P_value;
        else
            [~,ind21] = vlookupFast(x.rsID, PAR.RS);
            PAR.P = cat( 2, PAR.P, x.P_value(ind21));
        end
    end
    save([INPUT_DIR, 'par.mat'], 'PAR', '-v7.3');
end
%%
datasets = {'STAGE00DATA'};
datasetInd=1
    dataset = datasets{datasetInd};
    INPUT_DIR = ['../results/asymmetry/meta_analysis/' dataset '/not_imputed/not_subsampled/'];
    PAR = struct;
    for part=1:31
        disp(['Partition:', num2str(part)]);
        inpStatsFile = [INPUT_DIR 'CCAPart' sprintf('%02d', part), '.csv.gz'];
        unzippedFile = gunzip(inpStatsFile);
        x = readtable(unzippedFile{1},FileType="text",Delimiter=',');
        delete(unzippedFile{1});
        if part == 1
            PAR.RS = x.rsID;
            try
                PAR.POS = x.position;
            catch
            end
            PAR.CHR = x.chromosome;
            PAR.P = x.P_value;
        else
            [~,ind21] = vlookupFast(x.rsID, PAR.RS);
            PAR.P = cat( 2, PAR.P, x.P_value(ind21));
        end
    end
    save([INPUT_DIR, 'par.mat'], 'PAR', '-v7.3');
