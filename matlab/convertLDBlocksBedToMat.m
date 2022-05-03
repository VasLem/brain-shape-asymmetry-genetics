LD_BLOCKS = readtable('../SAMPLE_DATA/ld_blocks_euro.bed', 'FileType','text');

LDblocks = struct;
LDblocks.CHR =cellfun(@(x)str2double(x(4:end)), LD_BLOCKS.chr);
LDblocks.RANGES = table2array(LD_BLOCKS(:, 2:3));
REF.LDblocks = LDblocks;
save('../SAMPLE_DATA/ldblocks_euro.mat','REF', '-mat')