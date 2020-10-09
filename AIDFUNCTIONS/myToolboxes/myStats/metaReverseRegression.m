function [beta, Qe, pQe, Qb, pQb, mse, stdbeta] = metaReverseRegression(beta,covb)
        nDB = size(beta,2);
        nPred = size(beta{1},1);
        b = [];
        W = [];
        SIGMA = zeros(nPred,nPred);
        counter = 1;
        for db=1:1:nDB
            nPred = size(beta{db},1);
            if nPred > 0
                b = [b; beta{db}];
                W = [W; eye(nPred)]; %#ok<*AGROW>
                SIGMA(counter:counter+nPred-1,counter:counter+nPred-1) = covb{db};
                counter = counter+nPred;
            end
        end
        %SIGMA = diag(diag(SIGMA));% setting off diagonal elements to zero.

        [beta,stdbeta,mse,S] = lscov(W,b,SIGMA);

        Res = b - W*beta;
        Qe = Res'/SIGMA*Res;
        pQe = chi2cdf(Qe,(nDB-1)*nPred, 'upper');
        Qb = beta'/S*beta;
        pQb = chi2cdf(Qb, nPred, 'upper');
end


% function out = metaReverseRegression(beta,covb,x,y)
% 
%     nDB = size(beta,2);
%     nPred = size(beta{1},1);
%     b = [];
%     W = [];
%     SIGMA = zeros(nPred*nDB,nPred*nDB);
%     counter = 1;
%     for db=1:1:nDB
%         b = [b; beta{db}];
%         W = [W; eye(nPred)];
%         SIGMA(counter:counter+nPred-1,counter:counter+nPred-1) = covb{db};
%         counter = counter+nPred;
%     end
%     
%     res = regstats(b,W,'linear');
%     out.pFast = res.fstat.pval;
%     
% %     
% %     SIGMAInv = inv(SIGMA);
% %     covB = ((W'*SIGMAInv*W)^-1);
% %     B = covB*W'*SIGMAInv*b;
% %     
% %     pval = zeros(1,nDB);
% %     for db=1:nDB
% %         %db=1;
% %         tmpx = x{db};
% %         tmpy = y{db};
% %         tmpx = [ones(size(tmpx,1),1) tmpx];
% %         yhat = tmpx*B;
% %         
% %         residuals = tmpy - yhat;
% %         nobs = length(tmpy);
% %         p = length(B);
% %         dfe = nobs-p;
% %         dft = nobs-1;
% %         ybar = mean(tmpy);
% %         sse = norm(residuals)^2;    % sum of squared errors
% %         ssr = norm(yhat - ybar)^2;  % regression sum of squares
% %         sst = norm(tmpy - ybar)^2;     % total sum of squares;
% %         
% %         
% %         d = struct;
% %         d.sse = sse;
% %         d.dfe = dfe;
% %         d.dfr = p-1;
% %         d.ssr = ssr;
% %         d.f = (d.ssr/d.dfr)/(d.sse/d.dfe);
% %         pval(db) = 1-fcdf(d.f,d.dfr, d.dfe);
% 
% %     end
% %     
% %     
% %     
% %     
% %     
% %    diagB = diag(covB);
% % 
% %    Z = B./(sqrt(diagB));
% %    
% %    Qe = (b-W*B)'*SIGMAInv*(b-W*B);
% %    dfe = (nDB-1)*(nPred);
% %    pQe = 1-chi2cdf(Qe,dfe);
% %    
% %    Qb = B'*covB*B; 
% %    pQb = 1-chi2cdf(Qb,nPred);
% %    
% %    % returning to output
% %    out.B = B;
% %    out.covB = covB;
% %    out.Qe = Qe;
% %    out.pQe = pQe;
% %    out.Qb = Qb;
% %    out.pQb = pQb;
% %    out.pFast = pfast(pval);
% %    
% 
% end