function [out] = dGP3GroupANOVAAnalysis(G,dGP,t)
          index = find(~isnan(G));
          G = G(index);
          dGP = dGP(index);
          % Three group ANOVA
              XG = cell(size(G));
              XG(find(G==-1)) = {'AA'};
              XG(find(G==1)) = {'BB'};
              XG(find(G==0)) = {'AB'};
              out.ANOVA.XG = XG;
              [out.ANOVA.FA,out.ANOVA.pFA,out.ANOVA.ppermFA,out.ANOVA.T] = myPermAnova(XG,dGP,t);
end


