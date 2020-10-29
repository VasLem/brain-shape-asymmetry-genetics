function out = ShapeModelPermTest(Data,RefScan,t)
         if nargin < 3, t= 0;end
         % Building the original model
         [n,nV] = size(Data);
         SM = shapePCA;
         SM.RefScan = RefScan;
         getAverage(SM,Data');
         getModel(SM,Data');
         Exp = SM.Explained;
         CS = cumsum(SM.Explained);
         out.SM = SM;
         out.CS = CS;
         out.Exp = Exp;
         if t == 0, return; end
         CSCount = false(length(CS),t);
         ExpCount = false(length(CS),t);
         tic;
         parfor i=1:t
             Datafor = Data;
             for l=1:1:nV
                 index = randperm(n);
                 Datafor(:,l) = Datafor(index,l);
             end
             SMfor = shapePCA;
             SMfor.RefScan = RefScan;
             getAverage(SMfor,Datafor');
             getModel(SMfor,Datafor');
             CSfor = cumsum(SMfor.Explained);
             CSCount(:,i) = CSfor>=CS;
             ExpCount(:,i) = SMfor.Explained >= Exp;
         end
         toc;
         out.permCS = sum(CSCount,2)/t;
         out.permExp = sum(ExpCount,2)/t;
end