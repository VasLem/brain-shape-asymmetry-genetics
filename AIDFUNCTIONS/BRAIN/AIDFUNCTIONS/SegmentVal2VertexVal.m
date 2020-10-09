function [VAL,OUT] = SegmentVal2VertexVal(rend,SegmentVal)
         VAL = nan*zeros(size(rend.ULABELS));
         counter = 0;
         for i=1:rend.UHI.nLC
             if rend.UMASK(i)==0, continue;end
             counter = counter+1;
             [l,c] = Ind2LC(rend.UHI,i);
             index = find(rend.ULABELS(l,:)==c);
             VAL(l,index) = SegmentVal(counter);
         end
         if nargout==1, return; end
         % filling in the nan values;
         OUT = VAL;
         for l=1:1:size(OUT,1)
             TMP = VAL(1:l,:);
             values = TMP(end,:);
             counter = 1;
             while ~isempty(find(isnan(values)))&&counter<l
                    index = find(isnan(values));
                    values(index) = TMP(end-counter,index);
                    counter = counter+1;
             end
             OUT(l,:) = values;
         end
end




% rend = Render{1};