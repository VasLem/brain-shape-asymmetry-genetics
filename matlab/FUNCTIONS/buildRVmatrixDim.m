function out=buildRVmatrixDim(data,type,nDim)
% This function builds up the RV matrix across the population faces.
% Input:
% data: m x n, m participants, n corresponds to the (x,y,z); (t); (x,y,z,t) coordinates of
% the LMs.
% Output:
% out: is the nbpts x nbpts RV matrix, where nbpts is the number of LMs
% data = resTS; nDim =4; 
% data = res; nDim =3;
    if nargin<3, type = 'cov';end
    nInd = size(data,1);
    nbpts = size(data,2)/nDim;     
    data = reshape(data',nDim,nbpts,nInd);
    data = permute(data,[3 1 2]); % m x nDim x nbpts
    RV = zeros(nbpts,nbpts);
    disp('BUILDING RV MATRIX');
    ppb = ParforProgressbar(nbpts);
    parfor i=1:nbpts
        % i = 1; j = 2;
        tmp = zeros(1,nbpts);
        V1 = double(squeeze(data(:,:,i))); % data: nsub * (3, 1) * nVec
        for j=i:nbpts
            V2 = double(squeeze(data(:,:,j))); %#ok<PFBNS>
            covMatrix = cov([V1,V2]);
            switch lower(type)
                case 'cov'
                    tmp(j) = Escoufier_coefficient_fast_nDim(covMatrix);
                case 'cor'
                    tmp(j) = Escoufier_coefficient_fast_nDim(corrcov(covMatrix));
            end
        end
        RV(i,:) = tmp;
        ppb.increment();
    end
    delete(ppb);
    out = triu(RV,1)'+ triu(RV);
end
