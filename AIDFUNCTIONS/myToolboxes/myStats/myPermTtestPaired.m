function [T,pT,ppermT] = myPermTtestPaired(Y,t)  
         [~,pT,~,stats] = ttest(Y(:,1),Y(:,2));
         T = stats.tstat;
         % permuting
         if t==0, ppermT = 'not permuted'; return;  end
         n = size(Y,1);
         Tnull = zeros(1,t);
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
             [~,~,~,statsfor] = ttest(Yfor(:,1),Yfor(:,2));
             Tnull(i) = statsfor.tstat;       
         end
         ppermT = sum(Tnull>=T)/t;
end
