function preprocessGenome(PHENO_IID, CHR, GENO_PATH, MEDIAN_IMPUTE, BLOCK_SIZE, RESULTS_DIR, OUT_DIR, THREADS, SCRATCH_ROOT)
load(populateCCAWorkspace(RESULTS_DIR, SCRATCH_ROOT, CHR)); %#ok<LOAD>
disp("Loading PLINK preprocessed data..")
tic;
obj = SNPLIB();
obj.nThreads = THREADS;
[snps, samples] = obj.importPLINKDATA(GENO_PATH);
geno = obj.UnpackGeno();
clear obj;
% Show missing values
f = figure();
histogram(100 * (sum(geno==255) / size(geno, 1) ));
xlabel('Missing Values (%)')
ylabel("SNPs #")
saveas(f,[CHR_DIR,'missing_histogram.svg'])
close(f)
toc;
% Remove indels.
disp("Removing SNPs containing InDels..");
tic;
indelFilter = strlength(snps.(4))==1 & strlength(snps.(5))==1;
genoPruned = geno(:, indelFilter);
snpsPruned = snps(indelFilter, :);
clear geno
clear snps

if MEDIAN_IMPUTE
    disp("Imputing missing using median..")
    tic;
    parfor i=1:size(genoPruned,2)
        s = genoPruned(:, i);
        m =s~=255;
        medians(i) = round(median(s(m)));
        s(~m) = medians(i);
        genoPruned(:, i) = s;
    end
    [cnt_unique, unique_a] = hist(medians,single(unique(medians))); %#ok<HIST>
    imputes = array2table(cnt_unique);
    imputes.Properties.VariableNames=strsplit(num2str(unique_a));
    writetable(imputes,IMPUTE_INFO_OUT,'Delimiter',' ');
    toc;
end
if isfile(PLINK_DATA_INFO_OUT)
    load(PLINK_DATA_INFO_OUT, 'intervals', 'genoIndex')
    genoPruned = genoPruned(genoIndex, :);
else
    % Align phenotype with genotype
    disp("Aligning genotype to phenotype..");
    tic;
    iid = samples.IID;
    genoId = str2double(iid);
    phenoId = str2double(PHENO_IID);
    assignmentMatrix = (phenoId == genoId');
    [phenoIndex, genoIndex] = find(assignmentMatrix);
    iid = iid(genoIndex);
    genoPruned = genoPruned(genoIndex, :);
    assert(all(phenoIndex'==1:length(phenoIndex)));
    clear assignmentMatrix
    toc;
    % Some checks
    %regarding the fact that positions need to be sorted
    assert(all(sort(snpsPruned.POS) == snpsPruned.POS));
    % and that all phenotype ids  correspond to genotype ones
    assert(all((1:length(phenoId)) == phenoIndex'));
    % Partition SNPs
    disp("Computing Intervals..")
    tic;
    intervals = getIntervals(snpsPruned);
    toc;
    %%
    disp("Computing sample sizes per SNP..")
    tic;
    sampleSizes = getSampleSizes(genoPruned, intervals);
    toc;
    
    disp("Saving computed data..")
    tic;
    save(PLINK_DATA_INFO_OUT, 'snpsPruned', 'intervals', 'iid', 'genoIndex', 'sampleSizes', '-v7.3');
    toc;
end

total = size(intervals, 1);
blocksN = ceil(total/BLOCK_SIZE);

%%
disp(['Splitting SNPs into ' num2str(blocksN) ' blocks..'])
tic;
if isfolder(SCRATCH_GENE_DIR), rmdir(SCRATCH_GENE_DIR, 's'); end
mkdir(SCRATCH_GENE_DIR)
for blockCnt=1: blocksN
    minInd = (1 + (blockCnt - 1) * BLOCK_SIZE);
    maxInd= min(total, (blockCnt * BLOCK_SIZE));
    blockIntervals = intervals(minInd:maxInd,:);
    genoBlock = genoPruned(:, blockIntervals(1,1):blockIntervals(end,2));
    save([OUT_DIR num2str(blockCnt) '.mat'], 'genoBlock', 'blockIntervals', '-v6');
end
end


function intervals = getIntervals(snps, geneSetMethod, namedArgs)
% geneSetMethod:
%  slidingWindow: overlapping set of SNPS are taken based on the provided windowSize
%  perSnp: multiallelic SNPS are considered in the same CCA
if nargin < 2
    geneSetMethod = 'perSNP';
end
if nargin < 3
    namedArgs.windowSize = 20000;
end
switch geneSetMethod
    case 'slidingWindow'
        snpCnt = 1;
        cnt = 1;
        intSlidWindow = zeros(size(snps, 1),2);
        while snpCnt<size(snps, 1)
            startInd = snpCnt;
            startPos = snps.POS(snpCnt);
            offset = find(snps.POS(snpCnt:end) > startPos + namedArgs.windowSize, 1,  'first');
            if size(offset)==0, offset = size(snps, 1) + 1; end
            endInd = snpCnt + offset - 1;
            snpCnt = endInd + 1;
            intSlidWindow(cnt, :) = [startInd, endInd];
            cnt = cnt +1;
        end
        intervals = intSlidWindow(1:cnt-1, 2);
    case 'perSNP'
        [uniquePos, ia, ~] = unique(snps.POS);
        intervals = zeros(length(uniquePos), 2);
        intervals(:, 1) = ia(1:end);
        intervals(:, 2) = [ia(2:end)-1;length(snps.POS)];
end
end

function sampleSizes = getSampleSizes(geno, intervals)
t = sum(geno~=255);
parfor i=1:size(intervals,1)
    sampleSizes(i) = min(t(intervals(i,1): intervals(i,2)))
end
end