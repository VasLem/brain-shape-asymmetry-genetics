function [GC,pGC,seGC] = geneticCorrelationLDSCOV(CHR,POS,V1,V2,LDblocks,L2)
         nV1 = size(V1,2);
         nV2 = size(V2,2);
         disp('IDENTIFYING LD BLOCKS');
         LDblockID = getLDblockID(CHR,POS,LDblocks);
         Blocks = setdiff(unique(LDblockID),0);
         GC = zeros(nV1,nV2);
         pGC = ones(nV1,nV2);
         seGC = ones(nV1,nV2);
         disp('COMPUTING PAIRWISE CORRELATIONS');
         [path,ID] = setupParForProgress(nV1);
         parfor i=1:nV1
             % i=1;
             vv1 = V1(:,i);
             forGC = zeros(1,nV2);
             forpGC = ones(1,nV2);
             forseGC = ones(1,nV2);
             for j=1:nV2
                 % j=10; 
                 vv2 = V2(:,j); %#ok<*PFBNS>
                 [forGC(j),forseGC(j),forpGC(j)] = myLDSCOVwithSE(vv1(:),vv2(:),L2,LDblockID,Blocks);
             end
             GC(i,:) = forGC;
             pGC(i,:) = forpGC;
             seGC(i,:) = forseGC;
             parfor_progress;
         end
         closeParForProgress(path,ID)
end

function [stat,se,p] = myLDSCOVwithSE(Z1,Z2,L2,LDblockID,Blocks)
       Z1 = Z1(:);
       Z2 = Z2(:);
       L2= L2(:);
       nSNP = length(L2);
       W = ones(nSNP,1);
       X = lscov([ones(nSNP,1) L2 Z1],Z2,W);
       X1 = X(end);
       X = lscov([ones(nSNP,1) L2 Z2],Z1,W);
       X2 = X(end);
       stat = (X1+X2)/2; 
       n = length(Blocks);
       val = zeros(1,n);
       parfor i=1:n
           index = find(~(LDblockID==Blocks(i)));
           X = lscov([ones(length(index),1) L2(index) Z1(index)],Z2(index),W(index));
           X1 = X(end);
           X = lscov([ones(length(index),1) L2(index) Z2(index)],Z1(index),W(index));
           X2 = X(end);
           val(i) = (X1+X2)/2;   
       end
       se = std(val);
       p = 2*normcdf(-1*abs(stat/se),0,1);
end