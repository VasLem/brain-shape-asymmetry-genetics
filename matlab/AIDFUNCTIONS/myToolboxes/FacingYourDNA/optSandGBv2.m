function [optX,optE] = optSandGBv2(TMatches,FMatches,wM,nr,nSNP,X0)
         %X0 = ones(nr,1); 
         
         options = optimset('fminsearch');
         options.Display = 'none';
         
         nF = size(FMatches,3);
         nT = size(TMatches,2);
         
         
         T1Matches = TMatches(1:nr,:);
         T2Matches = TMatches(nr+1:end,:);
         F1Matches = FMatches(1:nr,:,:);
         F2Matches = FMatches(nr+1:end,:,:);
         
         wM1 = wM(1:nr,:);
         wM2 = wM(nr+1:end,:);
         
         T2Matches = nansum(wM2.*T2Matches,1);
         F2Matches = squeeze(nansum(repmat(wM2,1,1,nF).*F2Matches,1));
         
         
         TN = nansum(wM2,1);
         FN = repmat(TN,1,1,nF);
         
         
         [optX,optE] = fminsearch(@(X) myError(X,T1Matches,F1Matches,wM1,T2Matches,F2Matches,TN,FN,nT,nF),X0,options);
         optX = abs(optX);
         
end

function error = myError(X,T1Matches,F1Matches,wM1,T2Matches,F2Matches,TN,FN,nT,nF)
         [R,EER] = scoreAndEvalTF(T1Matches,F1Matches,wM1,T2Matches,F2Matches,TN,FN,nT,nF,abs(X));
         error = mean(R)+EER;
         
end

function [R,EER,ST,SF] = scoreAndEvalTF(T1Matches,F1Matches,wM1,T2Matches,F2Matches,TN,FN,nT,nF,W)   
         tmpW = repmat(W,1,nT).*wM1;  
         ST = ((nansum(tmpW.*T1Matches,1)+T2Matches)./(nansum(tmpW,1)+TN))';
         tmpW = repmat(tmpW,1,1,nF);
         SF = (squeeze(nansum(tmpW.*F1Matches))+F2Matches)./(squeeze(nansum(tmpW,1)+FN));
         [R,EER] = evalScores(ST,SF,nT,nF); 
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

