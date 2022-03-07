%runCCAonJobs

% Change the following
PHENO = 'pheno.mat';
JOBS_DIR = 'jobs';
RESULTS_DIR = 'gwas';
addpath(genpath('AIDFUNCTIONS'));
RESULTS_FILES = runCCA(PHENO, JOBS_DIR, RESULTS_DIR);