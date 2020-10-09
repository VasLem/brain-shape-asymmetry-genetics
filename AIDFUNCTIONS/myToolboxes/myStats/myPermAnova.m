function [FA,pFA,ppermFA,T] = myPermAnova(XG,Y,t)  
         [~,T] = anova1(Y,XG,'off');
         FA = T{2,5};
         pFA  = T{2,6};
         % permuting
         if t==0, ppermFA = nan; return;  end
         n = length(XG);
         FAnull = zeros(1,t);
         parfor i=1:t
             ind = randperm(n);
             Yfor = Y;
             Yfor = Yfor(ind,:);
             [~,Tfor] = anova1(Yfor,XG,'off');
             FAnull(i) = Tfor{2,5};       
         end
         ppermFA = sum(FAnull>=FA)/t;
end
