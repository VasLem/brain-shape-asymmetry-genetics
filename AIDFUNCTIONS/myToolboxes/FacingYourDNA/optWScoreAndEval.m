function [R,EER,S] = optWScoreAndEval(Matches,wM,W,n,g)
         S = getScores(Matches,wM,W,n);
         [R,EER] = evalScores(S,n,g); 
end

function Score = getScores(Matches,wM,W,n)
         Score = squeeze(nansum(repmat(W,1,n,n).*wM.*Matches,1));
         ScoreW = squeeze(nansum(repmat(W,1,n,n).*wM,1));
         Score = Score./ScoreW;
end

function [R,EER] = evalScores(Score,n,g)
         % Rank analysis
         R = mean((sum(Score<repmat(diag(Score),1,n),2)+1)/n);
         % EER analysis
         if nargout < 2, return; end
         [~,order] = sort(Score(:));
         g = g(order);
         true_neg = g == 0;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
         true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
         y = tpf(1:end-1)+dy./2;
         x = fpf(1:end-1)+dx./2;
         yn = 1-y;d = abs(x-yn);
         [~,indmin] = min(d);
         EER = ((x(indmin)+yn(indmin))/2);
end