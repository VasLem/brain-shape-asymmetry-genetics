function [out,index] = sortPeaks(PEAKS)
         c = unique(PEAKS.CHR);
         index = [];
         for i=1:length(c)
             indc = find(PEAKS.CHR==c(i));
             [~,sortind] = sort(PEAKS.POS(indc),'ascend');
             index = [index; indc(sortind)]; %#ok<*AGROW>
         end
         out = myReduceGWAS(PEAKS,index);
end