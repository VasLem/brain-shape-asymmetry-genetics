function tmppath=populateCCAWorkspace(RESULTS_DIR, SCRATCH_DIR, CHR)
tmppath=tempname;
SCRATCH_GENE_DIR = [SCRATCH_DIR 'jobs/' num2str(CHR) '/']; %#ok<NASGU> 
SCRATCH_OUTPUT_GENE_DIR = [SCRATCH_DIR 'jobs_output/' num2str(CHR) '/']; %#ok<NASGU> 
CHR_DIR = [RESULTS_DIR 'chr' num2str(CHR) '/'];
IMPUTE_INFO_OUT = [CHR_DIR 'imputed_freq.csv']; %#ok<NASGU>
PLINK_DATA_INFO_OUT = [CHR_DIR 'plink_data_info.mat']; %#ok<NASGU> 
META_INT_GENO_OUT = [CHR_DIR 'intervals_info.mat']; %#ok<NASGU> 
WITH_PART_CCA_OUT = [CHR_DIR 'withPartCCA.mat'];
WITH_PART_CCA_PROC = ~isfile(WITH_PART_CCA_OUT); %#ok<NASGU> 
SAMPLE_SIZES_OUT = [CHR_DIR 'sampleSizes.mat']; %#ok<NASGU>
if ~isfolder(CHR_DIR), mkdir(CHR_DIR); end
save(tmppath);
end