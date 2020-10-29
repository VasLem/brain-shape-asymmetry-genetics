function [LocalE,LocalS] = getLocalStatistic(M,R,B,type)
         [n,nB] = size(B); 
         P = B-R;
         A = repmat(mean(B),n,1);
         SST = sum((B-A).^2);
         SSR = sum((P-A).^2);
         switch lower(type)
             case 'value'
                 LocalE = M(2,:);
                 LocalS = SSR./SST;
              case 'shape'
                 LocalE = zeros(4,nB/3);
                 LocalE(1:3,:) = reshape(M(2,:),3,nB/3);
                 LocalE(4,:) = sqrt(sum(LocalE(1:3,:).^2));
                 SSR = reshape(SSR,3,nB/3);
                 SSR = sum(SSR);
                 SST = reshape(SST,3,nB/3);
                 SST = sum(SST);
                 LocalS = SSR./SST;
             otherwise
                 error('unknown type');
         end        

end
