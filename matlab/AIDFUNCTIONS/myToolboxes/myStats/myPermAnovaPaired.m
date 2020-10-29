function [FA,pFA,ppermFA,T] = myPermAnovaPaired(XG,Y,t)  
         [~,T] = anova1(Y,XG,'on');
         FA = T{2,5};
         pFA  = T{2,6};
         % permuting
         if t==0, ppermFA = 'not permuted'; return;  end
         n = size(Y,1);
         FAnull = zeros(1,t);
         Y1 = Y(:,1);
         Y2 = Y(:,2);
         parfor i=1:t
             r = randi(2,n,1);
             Y1perm = zeros(size(Y1));
             Y2perm = zeros(size(Y2));
             for k=1:1:n
                 switch r(k)
                     case 1
                        Y1perm(k,:) = Y1(k,:);
                        Y2perm(k,:) = Y2(k,:);
                      case 2
                         Y1perm(k,:) = Y2(k,:);
                         Y2perm(k,:) = Y1(k,:);
                 end
             end
             Yfor = [Y1perm Y2perm];
             [~,Tfor] = anova1(Yfor,XG,'off');
             FAnull(i) = Tfor{2,5};       
         end
         ppermFA = sum(FAnull>=FA)/t;
end
