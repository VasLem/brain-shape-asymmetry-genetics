function [out] = SubspaceDistance(X1,X2,center,AM,percvar,t)
    % This function computed the fluctuating direction Differences based on
    % unbalanced (unpaired) input arguments X1 & X2;
    if nargin < 4, t = 0; end
    if nargin < 2, center = true; end
    nX1 = size(X1,1);
    nX2 = size(X2,1);
    if percvar < 100
       disp('Computing Shape Model'); 
       SM = shapePCA;
       SM.RefScan = AM;
       getAverage(SM,X1');
       getModel(SM,X1');
       stripPercVar(SM,percvar);
       ncomp1 = SM.nrEV;
       SM = shapePCA;
       SM.RefScan = AM;
       getAverage(SM,X2');
       getModel(SM,X2');
       stripPercVar(SM,percvar);
       ncomp2 = SM.nrEV;
       ncomp = max(ncomp1,ncomp2);
    else
        ncomp = +inf;
    end
    %n = 50;
    AvgX1 = mean(X1);
    AvgX2 = mean(X2);
    if center
       X1 = X1-repmat(AvgX1,nX1,1);
       X2 = X2-repmat(AvgX2,nX2,1);
    end
    A = my_subspacea(X1',X2',ncomp);
    Dp = ProjectionMetric(A);
    Dpr = ProcrustesMetric(A);
    Dm = MitteroeckerMetric(A);
    out.A = A;
    out.Dprojection = Dp;
    out.Dprocrustes = Dpr;
    out.Dmitteroecker = Dm;
    DpCount = false(1,t);
    DprCount = false(1,t);
    DmCount = false(1,t);
    tic
    nT = nX1+nX2;
    X = [X1; X2];
    disp('Permuting');
    parfor i=1:t
        ind = randperm(nT);
        X1for = X(ind(1:nX1),:); %#ok<*PFBNS>
        X2for = X(ind(nX1+1:end),:);
        Afor = my_subspacea(X1for',X2for',ncomp);
        Dpfor = ProjectionMetric(Afor);
        Dprfor = ProcrustesMetric(Afor);
        Dmfor = MitteroeckerMetric(Afor);
        DpCount(i) = Dpfor>=Dp;
        DprCount(i) = Dprfor>=Dpr;
        DmCount(i) = Dmfor>=Dm;
    end
    disp('Done');
    toc;
    out.pDp = (sum(DpCount)/t);
    out.pDpr = (sum(DprCount)/t);
    out.pDm = (sum(DmCount)/t);
    out.ncomp = ncomp;
end

% function [out] = SubspaceDistance(X1,X2,center,AM,percvar,t)
%     % This function computed the fluctuating direction Differences based on
%     % unbalanced (unpaired) input arguments X1 & X2;
%     if nargin < 4, t = 0; end
%     if nargin < 2, center = true; end
%     nX1 = size(X1,1);
%     nX2 = size(X2,1);
%     if percvar < 100
%        disp('Computing Shape Model'); 
%        SM = shapePCA;
%        SM.RefScan = AM;
%        getAverage(SM,[X1',X2']);
%        getModel(SM,[X1',X2']);
%        stripPercVar(SM,percvar);
%        ncomp = SM.nrEV;
%     else
%         ncomp = +inf;
%     end
%     %n = 50;
%     AvgX1 = mean(X1);
%     AvgX2 = mean(X2);
%     if center
%        X1 = X1-repmat(AvgX1,nX1,1);
%        X2 = X2-repmat(AvgX2,nX2,1);
%     end
%     A = my_subspacea(X1',X2',ncomp);
%     Dp = ProjectionMetric(A);
%     Dpr = ProcrustesMetric(A);
%     Dm = MitteroeckerMetric(A);
%     out.A = A;
%     out.Dprojection = Dp;
%     out.Dprocrustes = Dpr;
%     out.Dmitteroecker = Dm;
%     DpCount = false(1,t);
%     DprCount = false(1,t);
%     DmCount = false(1,t);
%     tic
%     nT = nX1+nX2;
%     X = [X1; X2];
%     disp('Permuting');
%     parfor i=1:t
%         ind = randperm(nT);
%         X1for = X(ind(1:nX1),:); %#ok<*PFBNS>
%         X2for = X(ind(nX1+1:end),:);
%         Afor = my_subspacea(X1for',X2for',ncomp);
%         Dpfor = ProjectionMetric(Afor);
%         Dprfor = ProcrustesMetric(Afor);
%         Dmfor = MitteroeckerMetric(Afor);
%         DpCount(i) = Dpfor>=Dp;
%         DprCount(i) = Dprfor>=Dpr;
%         DmCount(i) = Dmfor>=Dm;
%     end
%     disp('Done');
%     toc;
%     out.pDp = (sum(DpCount)/t);
%     out.pDpr = (sum(DprCount)/t);
%     out.pDm = (sum(DmCount)/t);
%     out.ncomp = ncomp;
% end