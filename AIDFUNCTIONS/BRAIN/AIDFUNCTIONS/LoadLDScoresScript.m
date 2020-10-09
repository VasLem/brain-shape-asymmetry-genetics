function LDSCORE = LoadLDScoresScript
    %% script to load LD scores
    datapath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/LDSTUFF/ldsc-master/eur_w_ld_chr/';
    cd(datapath);
    bf = '.l2.ldscore';
    nCHR = 22;
    LDSCORE.CHR = [];
    LDSCORE.POS = [];
    LDSCORE.RS = {};
    LDSCORE.MAF = [];
    LDSCORE.L2 = [];
    for c=1:nCHR
        % c=22
        disp(num2str(c));
        file = [num2str(c) bf]; 
        T = readtable(file,'FileType','text');
        LDSCORE.CHR = [LDSCORE.CHR; uint8(T.CHR)];
        LDSCORE.POS = [LDSCORE.POS;uint32(T.BP)];
        LDSCORE.RS = [LDSCORE.RS;T.SNP(:)];
        LDSCORE.MAF = [LDSCORE.MAF;uint8(100*T.MAF)];
        %LDSCORE.L2 = [LDSCORE.L2;uint32(T.L2*10000)];
        LDSCORE.L2 = [LDSCORE.L2;T.L2];
    end
    %save('LDSCORE','LDSCORE','-v7.3');
end