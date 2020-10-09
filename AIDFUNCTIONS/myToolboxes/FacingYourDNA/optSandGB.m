function [optX,optE] = optSandGB(TMatches,FMatches,wM,nr,nSNP)
         X0 = ones(nr,1); 
         S0 = ones(nSNP,1);
         options = optimset('fminsearch');
         options.Display = 'none';
         [optX,optE] = fminsearch(@(X) myError(X,TMatches,FMatches,wM,S0),X0,options);
         optX = abs(optX);
         %[optX,optE] = fmincon(@(X) myError(X,TMatches,FMatches,wM,S0),X0,[],[],[],[],zeros(nr,1),[],[]);
end

function error = myError(X,TMatches,FMatches,wM,S0)
         [R,EER] = scoreAndEvalTF(TMatches,FMatches,wM,[abs(X);S0]);
         %[R,EER] = scoreAndEvalTF(TMatches,FMatches,wM,abs(X));
         error = mean(R)+EER;
         %error = EER;
         %error = mean(R);
end

function [R,EER,ST,SF] = scoreAndEvalTF(TMatches,FMatches,wM,W)
         nF = size(FMatches,3);
         nT = size(TMatches,2);
         tmpW = repmat(W,1,nT).*wM;
         ST = getScores(TMatches,tmpW)';
         SF = getScores(FMatches,repmat(tmpW,1,1,nF));
         [R,EER] = evalScores(ST,SF,nT,nF); 
         %R = 0;EER = 0;
end
function Score = getScores(Matches,tmpW)
         Score = squeeze(nansum(tmpW.*Matches,1)./nansum(tmpW,1));
end


function [R,EER] = evalScores(ST,SF,nT,nF)
         % Rank analysis
         R = (sum(SF<repmat(ST,1,nF),2)+1)/(nF+1);
         % EER analysis
         %if nargout < 2, return; end
         g = [ones(1,nT), zeros(1,nT*nF)];
         [~,order] = sort([ST;SF(:)]);
         g = g(order);
         true_neg = g == 0;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
         true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
         y = tpf(1:end-1)+dy./2;
         x = fpf(1:end-1)+dx./2;
         yn = 1-y;d = abs(x-yn);
         [~,indmin] = min(d);
         EER = ((x(indmin)+yn(indmin))/2);
end