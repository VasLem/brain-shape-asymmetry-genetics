function [GC,pGC,seGC] = geneticCorrelationLDblocksWithSE(CHR,POS,V1,V2,LDblocks,type)
         if nargin<6, type = 'mean';end
         nV1 = size(V1,2);
         nV2 = size(V2,2);
         disp('IDENTIFYING LD BLOCKS');
         LDblockID = getLDblockID(CHR,POS,LDblocks);
         Blocks = setdiff(unique(LDblockID),0);
         nBlocks = length(Blocks);
         v1 = zeros(nBlocks,nV1);
         v2 = zeros(nBlocks,nV2);
         %l2 = zeros(nBlocks,1);
         disp('DIVIDING INTO LD BLOCKS');
         for i=1:1:nBlocks
             switch type
                 case 'mean'
                    v1(i,:) = nanmean(V1(LDblockID==Blocks(i),:),1);
                    v2(i,:) = nanmean(V2(LDblockID==Blocks(i),:),1);
                 case 'median'
                    v1(i,:) = nanmedian(V1(LDblockID==Blocks(i),:),1);
                    v2(i,:) = nanmedian(V2(LDblockID==Blocks(i),:),1);
                 case 'max'
                    v1(i,:) = nanmax(V1(LDblockID==Blocks(i),:),1);
                    v2(i,:) = nanmax(V2(LDblockID==Blocks(i),:),1);
                 otherwise
                    v1(i,:) = nanmean(V1(LDblockID==Blocks(i),:),1);
                    v2(i,:) = nanmean(V2(LDblockID==Blocks(i),:),1);
             end
             %l2(i) = mean(L2(LDblockID==Blocks(i)));
         end
         GC = zeros(nV1,nV2);
         pGC = ones(nV1,nV2);
         seGC = ones(nV1,nV2);
         disp('COMPUTING PAIRWISE CORRELATIONS');
         [path,ID] = setupParForProgress(nV1);
         parfor i=1:nV1
             % i=1;
             vv1 = v1(:,i);
             forGC = zeros(1,nV2);
             forpGC = ones(1,nV2);
             forseGC = ones(1,nV2);
             for j=1:nV2
                 % j=200; 
                 vv2 = v2(:,j); %#ok<*PFBNS>
                 [forGC(j),forseGC(j),forpGC(j)] = getBootSTATS(vv1(:),vv2(:));
             end
             GC(i,:) = forGC;
             pGC(i,:) = forpGC;
             seGC(i,:) = forseGC;
             parfor_progress;
         end
         closeParForProgress(path,ID)
end
function [stat,se,p] = getBootSTATS(vv1,vv2)
         vv1 = double(vv1);
         vv2 = double(vv2);

         stat = corr(vv1(:),vv2(:),'type','Spearman');
         %stat = corr(vv1(:),vv2(:),'type','Pearson');
         n = length(vv1);
         val = zeros(1,100);
         for i=1:100
             index = randsample(n,n,true);
             val(i)= corr(vv1(index),vv2(index),'type','Spearman');
             %val(i)= corr(vv1(index),vv2(index),'type','Pearson');
         end
         se = std(val);
         %p = 2*normcdf(-1*abs(stat/se),0,1);
         p = normcdf(-1*(stat/se),0,1);      
end
% function [stat,se,p] = getJacknifeSTATS(vv1,vv2)
%          stat = corr(vv1(:),vv2(:),'type','Spearman');
%          n = length(vv1);
%          K = 20;
%          c = cvpartition(n,'KFold',K); 
%          val = zeros(1,K);
%          for i=1:K
%              index = training(c,i);
%              val(i)= corr(vv1(index),vv2(index),'type','Spearman');
%          end
%          se = std(val);
%          %p = 2*normcdf(-1*abs(stat/se),0,1);
%          p = normcdf(-1*(stat/se),0,1);
% end

