function [stat,statType,statP] = performStat(vIn,vripIn,var,t)
        stat = zeros(1,length(var.Booting));
        statType = cell(1,length(var.Booting));
        statP = zeros(1,length(var.Booting));
        for s=1:1:length(var.Booting)
           v = vIn(:,s);
           vrip = vripIn(:,s);
           index = find(~isnan(v));
           switch var.Info{s}.Type
               case 'Categorical'
                   statType{s} = 'F';
                   XG = cell(size(v(index)));
                   for l=1:length(var.Info{s}.El)
                       XG((v(index)==var.Info{s}.El(l))) = {num2str(l)};
                   end
                   [stat(s),~,statP(s)] = myPermAnova(XG,vrip(index),t);
               case 'Continous'
                   statType{s} = 'Cor';
                   [stat(s),statP(s)] = permCorr(v(index),vrip(index),t);
           end
        end  
end

        