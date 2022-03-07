#remotes::install_github("privefl/bigsnpr")
library(bigsnpr)
library(bigassertr)
system_verbose <- function(..., verbose) {
  system(..., ignore.stdout = !verbose, ignore.stderr = !verbose)
}

snp_beagleImpute <- function (beagle.path, plink.path, bedfile.in, bedfile.out = NULL, 
                              memory.max = 3, ncores = 1, extra.options = "", plink.options = "", 
                              verbose = TRUE) 
{
  prefix.in <- sub_bed(bedfile.in)
  if (is.null(bedfile.out)) 
    bedfile.out <- paste0(prefix.in, "_impute.bed")
  assert_noexist(bedfile.out)
  prefix.out <- sub_bed(bedfile.out)
  tmpfile1 <- tempfile()
  tmpfile2 <- tempfile()
  system_verbose(paste(plink.path, "--bfile", prefix.in, "--recode vcf bgz", 
                       "--out", tmpfile1, "--threads", ncores, plink.options), 
                 verbose = verbose)
  vcf1 <- paste0(tmpfile1, ".vcf.gz")
  on.exit(file.remove(vcf1), add = TRUE, after = FALSE)
  system_verbose(paste("java ", sprintf("-Xmx%dg", round(memory.max)), 
                       "-jar", beagle.path, paste0("gt=", vcf1), paste0("out=", 
                                                                        tmpfile2), paste0("nthreads=", ncores), extra.options), 
                 verbose = verbose)
  vcf2 <- paste0(tmpfile2, ".vcf.gz")
  on.exit(file.remove(vcf2), add = TRUE, after = FALSE)
  system_verbose(paste(plink.path, "--vcf", vcf2, "--out", 
                       prefix.out, plink.options), verbose = verbose)
  bedfile.out
}

for (DATASET in 1:2)
{
  if (DATASET==1){
    DATASET_ID='STAGE00DATA'
    BIOBANK_ID = 'UKBIOBANK'
    GENOME_ID = 'sel19908'
  }else{
    DATASET_ID = 'BATCH2_2021_DATA'
    BIOBANK_ID = 'MY_UKBIOBANK'
    GENOME_ID =  'sel16875_rmrel'
  }
  for (CHR in 1:22){
    print(sprintf('IMPUTING DATASET %d, CHROMOSOME %d', DATASET, CHR))
    input.file <- sprintf('../SAMPLE_DATA/IMAGEN/BRAIN/%s/GENOTYPES/PLINK/ukb_img_maf0.01_geno0.5_hwe1e-6_%s_chr%d.bed',BIOBANK_ID,GENOME_ID,CHR)
    output.dir <- sprintf('../SAMPLE_DATA/IMAGEN/BRAIN/IMPUTED_GENOTYPES/%s',DATASET_ID)
    dir.create(output.dir,showWarnings = FALSE,recursive=T)
    output.file <- sprintf('%s/ukb_img_maf0.01_geno0.5_hwe1e-6_%s_beagle_chr%d.bed',output.dir, GENOME_ID,CHR)
    
    snp_beagleImpute(beagle.path='../binaries/beagle.08Feb22.fa4.jar',
                              plink.path='../binaries/plink/plink', 
                              bedfile.in=input.file,
                              bedfile.out=output.file,
                              ncores=8,
                              memory.max=16, extra.options = sprintf('window=20 map=../SAMPLE_DATA/genetic_map/plink.chr%d.GRCh37.map ref=../SAMPLE_DATA/1000GenomeProject/euro/chr%d.1kg.phase3.v5a.b37.euro.bref3',CHR, CHR)
    )
  }
}




