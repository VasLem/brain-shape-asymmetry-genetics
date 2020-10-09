function out = ShapeModelParallAnalysis(Data,RefScan,t)
         if nargin < 3, t= 0;end
         % Building the original model
         [n,nV] = size(Data);
         mu = mean(Data);
         sigma = std(Data); 
         SM = shapePCA;
         SM.RefScan = clone(RefScan);
         getAverage(SM,Data');
         getModel(SM,Data');
         Exp = SM.Explained;
         CS = cumsum(SM.Explained);
         EigVal = SM.EigVal;
         out.SM = SM;
         out.CS = CS;
         out.Exp = Exp;
         if t == 0, return; end
         CSCount = false(length(CS),t);
         ExpCount = false(length(CS),t);
         EigValCount = false(length(CS),t);
         tic;
         parfor i=1:t
             Datafor = repmat(mu,n,1) + repmat(sigma,n,1).*randn(n,nV);
             SMfor = shapePCA;
             SMfor.RefScan = clone(RefScan);
             getAverage(SMfor,Datafor');
             getModel(SMfor,Datafor');
             CSfor = cumsum(SMfor.Explained);
             CSCount(:,i) = CSfor>=CS;
             ExpCount(:,i) = SMfor.Explained >= Exp;
             EigValCount(:,i) = SMfor.EigVal >= EigVal;
         end
         toc;
         out.paCS = sum(CSCount,2)/t;
         out.nCS = find(out.paCS<=0.05,1,'last');
         out.paExp = sum(ExpCount,2)/t;
         out.nExp = find(out.paExp<=0.05,1,'last');
         out.paEigVal = sum(EigValCount,2)/t;
         out.nEigVal = find(out.paEigVal<=0.05,1,'last');
         % reducing the model
         reduceNrPc(SM,out.nEigVal);
         out.TotalExplained = CS(out.nEigVal);
end