function [FA,pFA,ppermFA,T] = myPermAnovav2(XG,Y,t)  
         [~,T] = anova1(Y,XG,'on');
         FA = T{2,5};
         pFA  = T{2,6};
         % permuting
         if t==0, ppermFA = 'not permuted'; return;  end
         n = size(Y,1);
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
