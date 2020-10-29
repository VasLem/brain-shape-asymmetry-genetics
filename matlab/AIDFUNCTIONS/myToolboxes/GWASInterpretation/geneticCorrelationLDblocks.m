function [GC,pGC,permGC] = geneticCorrelationLDblocks(CHR,POS,V1,V2,LDblocks,type,t)
         if nargin<7, t=0; end
         if nargin<6, type = 'mean';end
         nV1 = size(V1,2);
         nV2 = size(V2,2);
         disp('IDENTIFYING LD BLOCKS');
         LDblockID = getLDblockID(CHR,POS,LDblocks);
         Blocks = setdiff(unique(LDblockID),0);
         nBlocks = length(Blocks);
         v1 = zeros(nBlocks,nV1);
         v2 = zeros(nBlocks,nV2);
         disp('DIVIDING INTO LD BLOCKS');
         for i=1:1:nBlocks
             switch type
                 case 'mean'
                    v1(i,:) = mean(V1(LDblockID==Blocks(i),:));
                    v2(i,:) = mean(V2(LDblockID==Blocks(i),:));
                 case 'median'
                    v1(i,:) = median(V1(LDblockID==Blocks(i),:));
                    v2(i,:) = median(V2(LDblockID==Blocks(i),:));
                 case 'max'
                    v1(i,:) = max(V1(LDblockID==Blocks(i),:));
                    v2(i,:) = max(V2(LDblockID==Blocks(i),:));
                 otherwise
                    v1(i,:) = mean(V1(LDblockID==Blocks(i),:));
                    v2(i,:) = mean(V2(LDblockID==Blocks(i),:));
             end         
         end
         GC = zeros(nV1,nV2);
         pGC = ones(nV1,nV2);
         permGC = ones(nV1,nV2);
         disp('COMPUTING PAIRWISE CORRELATIONS');
         [path,ID] = setupParForProgress(nV1);
         parfor i=1:nV1
             vv1 = v1(:,i);
             forGC = zeros(1,nV2);
             forpGC = ones(1,nV2);
             forpermGC = ones(1,nV2);
             for j=1:nV2
                 vv2 = v2(:,j); %#ok<*PFBNS>
                 [forGC(j),forpGC(j)] = corr(vv1(:),vv2(:),'type','Spearman');
                 if t>0
                    forpermGC(j) = permTest(vv1(:),vv2(:),forGC(j),t);
                 end
             end
             GC(i,:) = forGC;
             pGC(i,:) = forpGC;
             permGC(i,:) = forpermGC;
             parfor_progress;
         end
         closeParForProgress(path,ID)
end
function permp = permTest(v1,v2,obsval,t)
         permval = zeros(1,t);
         %disp('PERMUTING');
         for i=1:t
             ind = randperm(length(v1));
             permval(i) = corr(v1(ind),v2,'type','Spearman');
         end
         permp = (sum(permval>=obsval)+1)/(t+1);
         %disp('DONE');
end