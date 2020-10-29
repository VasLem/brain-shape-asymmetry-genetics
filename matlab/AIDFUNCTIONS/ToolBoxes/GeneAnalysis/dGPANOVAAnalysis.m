function [out,fig] = dGPANOVAAnalysis(G,dGP,t,name)
%           G = GM(:,26);
%           dGP = dGPM(:,26);
%           t = 10000;
          index = find(~isnan(G));
          G = G(index);
          dGP = dGP(index);
          % Three group ANOVA
              disp('grouped anova');
              XG = cell(size(G));
              XG(find(G==-1)) = {'AA'};
              XG(find(G==1)) = {'BB'};
              XG(find(G==0)) = {'AB'};
              out.ANOVA.XG = XG;
              [out.ANOVA.FA,out.ANOVA.pFA,out.ANOVA.ppermFA,out.ANOVA.T] = myPermAnova(XG,dGP,t);
          % Two group ANOVA 1
              disp('Pair wise anova 1');
              index = find(G==1);
              ind = (1:size(G,1));
              ind2 = setdiff(ind,index);
              G2 = G(ind2);
              XG2 = cell(size(G2));
              XG2(find(G2==-1)) = {'AA'};
              XG2(find(G2==0)) = {'AB'};
              out.ANOVAPW(1).XG = XG2;
              out.ANOVAPW(1).Grouping = 'AA AB';
              [out.ANOVAPW(1).FA,out.ANOVAPW(1).pFA,out.ANOVAPW(1).ppermFA,out.ANOVAPW(1).T] = myPermAnova(XG2,dGP(ind2),t);
          % Two group ANOVA 2
              disp('Pair wise anova 2');
              index = find(G==-1);
              ind = (1:size(G,1));
              ind2 = setdiff(ind,index);
              G2 = G(ind2);
              XG2 = cell(size(G2));
              XG2(find(G2==1)) = {'BB'};
              XG2(find(G2==0)) = {'AB'};
              out.ANOVAPW(2).XG = XG2;
              out.ANOVAPW(2).Grouping = 'BB AB';
              [out.ANOVAPW(2).FA,out.ANOVAPW(2).pFA,out.ANOVAPW(2).ppermFA,out.ANOVAPW(2).T] = myPermAnova(XG2,dGP(ind2),t);
          % Two group ANOVA 1
              disp('Pair wise anova 3');
              index = find(G==0);
              ind = (1:size(G,1));
              ind2 = setdiff(ind,index);
              G2 = G(ind2);
              XG2 = cell(size(G2));
              XG2(find(G2==-1)) = {'AA'};
              XG2(find(G2==1)) = {'BB'};
              out.ANOVAPW(3).XG = XG2;
              out.ANOVAPW(3).Grouping = 'AA BB';
              [out.ANOVAPW(3).FA,out.ANOVAPW(3).pFA,out.ANOVAPW(3).ppermFA,out.ANOVAPW(3).T] = myPermAnova(XG2,dGP(ind2),t);
          % Generating Boxplot
              fig = [];
%               indAA = find(strcmp(XG,'AA'));
%               indBB = find(strcmp(XG,'BB'));
%               indAB = find(strcmp(XG,'AB'));
%               ind = [indAA(:); indAB(:); indBB(:)];          
%               %fig = figure;boxplot(dGP(ind),XG(ind));title(name); 
%               str = ['pAll:' num2str(out.ANOVA.ppermFA,'%1.3f') ...
%                      ' pAA2AB:' num2str(out.ANOVAPW(1).ppermFA,'%0.3f') ...
%                      ' pBB2AB:' num2str(out.ANOVAPW(2).ppermFA,'%0.3f') ...
%                      ' pAA2BB:' num2str(out.ANOVAPW(3).ppermFA,'%0.3f')];
%               %xlabel(str);
%               %drawnow;
end


