LD_BLOCKS = readtable('../SAMPLE_DATA/ld_blocks_euro.bed', 'FileType','text');
LDBlocks = struct;
LDblocks.CHR = table2array(LD_BLOCKS.chr);
LDBlocks.RANGE = table2array(LDBlocks(:, 2:3));
REF.LDBlocks = LDBlocks;
save('ldblocks_euro.mat','REF', '-mat')