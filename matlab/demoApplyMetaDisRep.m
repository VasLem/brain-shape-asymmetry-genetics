%ENV TO SET:
%%%%%%%%%%%%
%DATA_ROOT: the directory of the DATA (../SAMPLE_DATA/)
%RESULTS_ROOT: the directory of the results (../results/)
%SCRATCH_ROOT: the directory  of the temporary results (../results/)
%THREADS: Number of threads to use (max set by local)
%CHROMOSOME: Chomosome to analyze (1:22)
%BLOCK_SIZE: Block Size for CCA (2000)
%IMPUTE_STRATEGY: Whether to perform imputation (zero,mean,median or beagle) or not (no). No imputation is quite slow. (no)
%SUBSAMPLED: Whether to use the subsampled phenotype, if no PHENO_PATH is provided. It modifies the saving and scratch directories (defaults to 0)
%PHENO_PATH: Whether to use a specific path for the mat file of the phenotype (defaults to the path where hierarchical clustering algorithm places it)
%%%%%%%%%%%%
close all;clear;
if ~isdeployed
    restoredefaultpath;
    addpath(genpath('.'));
    addpath(genpath('AIDFUNCTIONS'));
    addpath(genpath('../SNPLIB/'));
end
%%
UPDATE_FIGS = 1;
NO_PARTITION_THRES = 5*10^-8; % European in LD score
DEFAULT_CHRS = 1:22;
DEFAULT_SUBSAMPLED = 0;
DEFAULT_MEDIAN_IMPUTE = 'mean';
DEFAULT_MODALITY = 'asymmetry'; %symmetry,asymmetry

SUBSAMPLED = getenv('SUBSAMPLED');
if(isempty(SUBSAMPLED))
    SUBSAMPLED = DEFAULT_SUBSAMPLED;
else
    if ~isnumeric(SUBSAMPLED)
        SUBSAMPLED=str2double(SUBSAMPLED);
    end
end
if SUBSAMPLED
    REDUCTION_ID='subsampled';
else
    REDUCTION_ID='not_subsampled';
end

IMPUTE_STRATEGY = getenv('IMPUTE_STRATEGY');
if(isempty(IMPUTE_STRATEGY))
    IMPUTE_STRATEGY = DEFAULT_MEDIAN_IMPUTE;
end
switch IMPUTE_STRATEGY
    case 'no'
        IMPUTE_ID = 'not_imputed';
    case 'zero'
        IMPUTE_ID = 'zero_imputed';
    case 'mean'
        IMPUTE_ID = 'mean_imputed';
    case 'median'
        IMPUTE_ID = 'median_imputed';
    case 'beagle'
        IMPUTE_ID = 'beagle_imputed';
    otherwise
        error("IMPUTE_STRATEGY not understood, available options: no zero mean median beagle")
end



DATA_DIR = getenv('DATA_ROOT');
if(isempty(DATA_DIR))
    DATA_DIR = '../SAMPLE_DATA/';
end

SCRATCH_ROOT = getenv('SCRATCH_ROOT');
if(isempty(SCRATCH_ROOT))
    SCRATCH_ROOT = '../results/';
end


THREADS = getenv('THREADS');
if(isempty(THREADS))
    THREADS = parcluster('local');
    THREADS = THREADS.NumWorkers;
else
    if ~isnumeric(THREADS)
        THREADS=str2double(THREADS);
    end
end

CHRS = getenv("CHROMOSOME");
if(isempty(CHRS))
    CHRS = DEFAULT_CHRS;
else
    if ~isnumeric(CHRS)
        CHRS=str2double( strsplit(CHRS,','));
    end
end

MODALITY = getenv("MODALITY");
if(isempty(MODALITY))
    MODALITY = DEFAULT_MODALITY;
end

disp(['Modality: ' MODALITY])
disp(['Number of threads: ', num2str(THREADS)])
disp(['Location of data: ', DATA_DIR])
disp(['Location of temporary data:', SCRATCH_ROOT])
disp(['Imputation: ', num2str(IMPUTE_STRATEGY)])


DATASET_ID = [IMPUTE_ID '/' REDUCTION_ID];
OUTPUT_DIR = ['../results/' MODALITY '/meta_analysis/joinedDatasets/' DATASET_ID '/'];
if ~isfolder(OUTPUT_DIR), mkdir(OUTPUT_DIR); end
%%
obj = SNPLIB();
obj.nThreads = THREADS;
disp("Reading SNPs info..")
for CHR=1:22
    GENO_PATH = [DATA_DIR 'IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_sel19908_chr' num2str(CHR)];
    [snps{CHR},~] = obj.importPLINKDATA(GENO_PATH);
end
disp("Joining datasets..")
snps = cat(1,snps{:});
for partition=1:31
    disp(['Partition ' num2str(partition)]);
    ftab = [];
    for datasetInd=1:2
        if datasetInd==1
            DATASET_NAME = 'STAGE00DATA';
        else
            DATASET_NAME = 'BATCH2_2021_DATA';
        end
        INPUT_DIR = ['../results/' MODALITY '/meta_analysis/' DATASET_NAME '/' DATASET_ID '/'];
        fname = gunzip(sprintf("%sCCAPart%02d.csv.gz",INPUT_DIR, partition));
        outname = sprintf("%sCCAPart%02d.csv",OUTPUT_DIR, partition);
        tab = readtable(fname{1});
        delete(fname{1});
        if isempty(ftab)
            ftab = tab;
        else
            ftab = innerjoin(ftab, tab, keys='rsID'); % overlapping SNPs only
            ftab = innerjoin(ftab, snps,  LeftKeys='rsID', RightKeys='RSID'); % get Positions
            pvalues = ftab(:,{'P_value_tab','P_value_ftab'});
            ftab.P_value = stouffer(table2array(pvalues));
            zeroMask = ftab.P_value == 0;
            ftab.P_value(zeroMask) = min(ftab.P_value(~zeroMask))/10;
            ftab.ChiScore = (ftab.ChiScore_tab + ftab.ChiScore_ftab) / 2;
            ftab.N = ftab.N_ftab + ftab.N_tab;
            ftab = ftab(:, {'rsID', 'POS', 'N','A2_ftab','A1_ftab', 'P_value','ChiScore', 'chromosome_ftab'});
            ftab.Properties.VariableNames = {'rsID','position', 'N','A2','A1','P_value','ChiScore','chromosome'};
            writetable(ftab,outname);
            gzip(outname);
            delete(outname);
        end

    end
end



