function [FACIALCOV,B] = learnFacialResidus(COV,SHAPE,K)
         if nargin<3, K = 10;end
         [nS,nC] = size(COV);
         [nV,~] = size(SHAPE);
         F = crossvalind('Kfold',nS,K);
         B = nan*zeros(nC,nV,K);
         FACIALCOV = nan*zeros(nS,nC);
         for k=1:1:K
             %k=1;
             Testing = find(F==k);
             Training = setdiff(1:nS,Testing);
             [~,~,~,~,betha] = plsregress(COV(Training,:),SHAPE(:,Training)');
             
             in = SHAPE(:,Training(1))';
             avgX = mean(COV(Training,:));
             
             
             X = avgX;
             M = betha;
             Y = in;
             
             
             B(:,:,k) = betha(2:end,:); 
             tmp = nan*zeros(length(Testing),nC);
             for c=1:1:nC
                 %c=1;
                 M = betha(1+c,:);
                 in = SHAPE(:,Testing)';
                 out = dot(in',repmat(M'/norm(M'),1,size(in,1)));
                 tmp(:,c) = out';
             end
             FACIALCOV(Testing,:) = tmp;
         end
end


