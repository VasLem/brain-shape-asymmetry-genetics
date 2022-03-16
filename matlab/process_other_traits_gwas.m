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
    RESULTS_DIR = ['../results/ldsc/' TRAIT_ID '/munged'];
    if ~isfolder(RESULTS_DIR), mkdir(RESULTS_DIR); end
    for partition=1:N_PARTITIONS
        ftab = table;
        ftab.SNP = struct.RS;
        ftab.N = struct.N;
        ftab.A1 = struct.A1;
        ftab.A2 = struct.A2;
        ftab.Z = struct.CHI(:, partition);
        fname =  [RESULTS_DIR '/' sprintf('par%02d.sumstats',partition)];
        writetable(ftab,fname,FileType='text', Delimiter='\t');
        gzip(fname);
        delete(fname);
    end
end