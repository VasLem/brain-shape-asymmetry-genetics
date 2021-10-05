addpath(genpath('/opt/SNPLIB/'));
chr = 1;

obj = SNPLIB();
obj.nThreads = 8;
[snps, samples] = obj.importPLINKDATA(['/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/PLINKSELECTEDFILTERED/ukb_img_maf0.01_geno0.5_hwe1e-6_sel19908_chr' num2str(chr)]);
geno = obj.UnpackGeno();