addpath(genpath('/opt/SNPLIB/'));
chr = 1;

obj = SNPLIB();
obj.nThreads = 8;
[snps, samples] = obj.importPLINKDATA(['../SAMPLE_DATA/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_sel19908_chr' num2str(chr)]);
geno = obj.UnpackGeno();
%%
af = obj.CalcAlleleFrequency();
%%
grm = obj.CalcGRMMatrix();