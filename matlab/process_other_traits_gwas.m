% load face and brain shape
face = load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/ANALYSIS/PEAKDETECTION/FACESUG.mat');
N_FACE_PARTITIONS = 63;
brain = load('/usr/local/micapollo01/IMAGEN_DATA/SHARED/pclaes4/UKB/DATA/ANALYSIS/PEAKDETECTION/BRAINSUG_SWH.mat');
N_BRAIN_PARTITIONS = 285;
%%

for c=1:2
    if c==1
        TRAIT_ID = 'face';
        N_PARTITIONS = 63;
        struct = face.FACESUG;
    else
        TRAIT_ID = 'brain_shape';
        N_PARTITIONS = 285;
        struct = brain.BRAINSUG;
    end
    RESULTS_DIR = ['../results/other_traits_gwas/' TRAIT_ID '/'];
    if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
    for partition=1:N_PARTITIONS
        ftab = table;
        ftab.rsID = struct.RS;
        ftab.N = struct.N;
        ftab.A1 = struct.A1;
        ftab.A2 = struct.A2;
        ftab.chromosome = struct.CHR;
        ftab.('P-value') = struct.P(:, partition);
        if c == 1
            ftab.ChiScore = struct.CHI(:, 1 + 2*(partition-1));
        else
            ftab.ChiScore = struct.CHI(:, partition);
        end
        fname =  [RESULTS_DIR '/' sprintf('CCAPart%02d.csv',partition)];
        writetable(ftab,fname);
        gzip(fname);
        delete(fname);
    end
end