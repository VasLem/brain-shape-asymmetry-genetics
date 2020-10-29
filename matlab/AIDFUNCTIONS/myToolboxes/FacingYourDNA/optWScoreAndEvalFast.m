function [R,EER,S] = optWScoreAndEvalFast(Matches,wM,W,n,nrW,Tind,Find)
         if nargin < 4
            n = size(Matches,2);
            nrW = size(W,2);
            ind1 = repmat((1:n)',1,nrW);ind1 = ind1(:);
            ind2 = repmat((1:nrW),n,1);ind2 = ind2(:);
            Aind = 1:(n*n*nrW);
            Tind = sub2ind([n n nrW],ind1,ind1,ind2);
            Find = setdiff(Aind,Tind);
         end
         S = getScores(Matches,wM,W,n,nrW);
         [R,EER] = evalScores(S,Tind,Find,n,nrW); 
end

function Score = getScores(Matches,wM,W,n,nrW)
         MS = repmat(wM.*Matches,1,1,1,nrW);
         WS = permute(repmat(W,1,1,n,n),[1 3 4 2]);
         Score = squeeze(nansum(WS.*MS,1));
         Score = Score./squeeze(nansum(WS.*repmat(wM,1,1,1,nrW),1));
end

function [R,EER] = evalScores(Score,Tind,Find,n,nrW)
         % Rank analysis   
         TS = reshape(Score(Tind),n,nrW);
         R = mean((squeeze(sum(Score<permute(repmat(TS,1,1,n),[1 3 2]),2))+1)/n,1);
         if nargout < 2, return; end
         FS = reshape(Score(Find),n*n-n,nrW);
         g = repmat([ones(1,n) zeros(1,n*(n-1))]',1,nrW);
         % EER analysis
         [~,order] = sort([TS;FS]);
         g = g(order);
         true_neg = g == 0;true_pos = g == 1;
         fpf = cumsum(true_neg,1)./(n*(n-1));tpf = cumsum(true_pos,1)./n;
         dx = diff(fpf);dy = diff(tpf);
         y = tpf(1:end-1,:)+dy./2;
         x = fpf(1:end-1,:)+dx./2;
         yn = 1-y;d = abs(x-yn);
         [~,indmin] = min(d);
         EER = ((x(indmin)+yn(indmin))/2);
end