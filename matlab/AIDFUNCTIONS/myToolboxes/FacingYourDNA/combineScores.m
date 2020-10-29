function out = combineScores(input,type)
         switch lower(type)
             case 'mean'
                 out = nanmean(input,1);
             case 'median'
                 out = nanmedian(input,1);
             case 'max'
                 out = nanmax(input,[],1);
             case 'min'
             case 'product'
         end
end