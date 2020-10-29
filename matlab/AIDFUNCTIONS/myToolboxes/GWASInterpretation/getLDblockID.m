function LDblockID = getLDblockID(chr,pos,LDblocks)
         nSNP = length(chr);
         if ~length(pos)==nSNP, error('CHR and POS should contain the same number of entries');end
         LDblockID = zeros(1,nSNP,'uint16');
         for i=1:LDblocks.n
             index = find(chr==LDblocks.CHR(i));
             if isempty(index), continue; end
             startrange = pos(index)>=LDblocks.RANGES(i,1);
             endrange = pos(index)<=LDblocks.RANGES(i,2);
             total = startrange+endrange;
             index2 = find(total==2);
             if isempty(index2), continue;end
             LDblockID(index(index2)) = i;
         end
end