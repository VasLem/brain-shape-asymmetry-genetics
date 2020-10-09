function out = optimizeSNPMatchWeights(testM,M,options)
         nSNP = length(testM);
         out = nan*zeros(1,nSNP);
         index = find(~isnan(testM));
         out(index) = 1; %#ok<*FNDSB>
         testM = -log(testM(index))';M = -log(M(index,:))';
         out(index) = getW(testM,M,options);
%          [bestW,avgW] = getW(testM,M,options);
%          outBest = out;outBest(index) = bestW;
%          outAvg = out;outAvg(index) = avgW;
end

function [bestW,avgW] = getW(testM,M,options)
         % initialize
         nW = length(testM);
         iter = 0;
         redind = ceil(options.Red*options.PopSize);
         % initialize
         Pop = rand(options.PopSize,nW);
         ranks = evalPop(Pop,testM,M);
         %disp(['Iter: ' num2str(0) ' Best: ' num2str(min(ranks)) ' Mean: ' num2str(mean(ranks))]);
         iterate = true;
         while iter<options.MaxIter&&iterate
            iter = iter + 1;
           % SORT
            [ranks,sortind] = sort(ranks); 
           % GENERATE
            Pop = nextPop(Pop,sortind,redind);
           % EVAL
            ranks(redind+1:end) = evalPop(Pop(redind+1:end,:),testM,M);
            %disp(['Iter: ' num2str(iter) ' Best: ' num2str(min(ranks)) ' Mean: ' num2str(mean(ranks))]);
            if length(find(ranks==0))>= options.minRank0, iterate = false; end
         end
         minrank = min(ranks);
         bestW = mean(Pop(ranks==minrank,:));
         %disp(num2str(length(find(ranks==minrank))));
         if nargout < 2, return; end
         [~,ind] = sort(ranks);
         scaling = zeros(size(ranks));
         scaling(ind) = 1./ ((1:length(ranks)).^0.5);
         scaling = options.PopSize * scaling ./ sum(scaling);
         avgW = sum(repmat(scaling,1,nW).*Pop)./sum(scaling);
end

function ranks = evalPop(Pop,testM,M)
         nrM = size(M,1);N = sum(Pop,2);PopSize = size(Pop,1);
         testS = (sum(Pop.*repmat(testM,PopSize,1),2)./N);
         S = squeeze(sum(repmat(Pop,1,1,nrM).*permute(repmat(M,1,1,PopSize),[3 2 1]),2))./repmat(N,1,nrM);
         ranks = sum(S<repmat(testS,1,nrM),2);
end

function Pop = nextPop(Pop,sortind,redind)
         Pop = Pop(sortind,:);PopSize = size(Pop,1);
         PopR = Pop(1:redind,:);
         RemainSize = PopSize - redind;
         h = 1.06*min(iqr(PopR)/1.349,std(PopR))./(redind^(1/5));
         PopR = PopR';c= size(PopR,2);
         randsel = bsxfun(@(x,y) x(randperm(y(1))),PopR',repmat(c,1,c)');
         Pop(redind+1:end,:) = randsel(1:RemainSize,:) + repmat(h,RemainSize,1).*randn(RemainSize,size(Pop,2));
         Pop(Pop<0) = 0;% exclude negative weigths
         %Pop(Pop>1) = 1;% upper boundary on weigths
end


% function [out,optRank] = optimizeSNPMatchWeights(testM,M,options)
%          nSNP = length(testM);
%          out = nan*zeros(1,nSNP);
%          index = find(~isnan(testM));
%          out(index) = 1; %#ok<*FNDSB>
%          testM = -log(testM(index));M = -log(M(index,:));
%          X0 = ones(size(testM));
%          nr = size(M,2)+1;
%          %X0 = rand(length(testM),1);
%          [optW,optRank] = getW(testM,M,nr,X0,options);
%          out(index) = optW;
% end


% function [optX,optRank] = getW(testM,M,nr,X0,options)
%          % Constrained optimization
%          %[optX,optRank] = fmincon(@(X) getRank(X,testM,M,nr),X0,[],[],[],[],zeros(size(testM)),2*ones(size(testM)),[],options);
%          
%          [optX,optRank] = ga(@(X) getRank(X,testM,M,nr),length(X0),[],[],[],[],zeros(size(testM)),2*ones(size(testM)),[],options); 
%          
%          
% end
% 
% function rank = getRank(X,testM,M,nr)        
%          N = sum(X);
%          testS = sum(X.*-log(testM))/N;
%          S = sum(repmat(X,1,size(M,2)).*-log(M))./N;
%          rank = (length(find(S<testS))+1)/nr;
% end
