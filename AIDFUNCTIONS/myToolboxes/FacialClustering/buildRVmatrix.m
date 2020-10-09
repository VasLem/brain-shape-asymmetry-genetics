function out=buildRVmatrix(data,type)

% This function builds up the RV matrix across the population faces.
% Input:
% data: m x n, m participants, n corresponds to the (x,y,z) coordinates of
% the LMs.
% Output:
% out: is the nbpts x nbpts RV matrix, where nbpts is the number of LMs
    if nargin<2, type = 'cov';end
    nbpts = size(data,2)/3;
    RV = zeros(nbpts,nbpts);
    disp('BUILDING RV MATRIX');
    [path,ID] = setupParForProgress(nbpts);
    parfor i=1:nbpts
        % i = 1;
        tmp = zeros(1,nbpts);
        for j=i:nbpts
            idx1 = repmat(((i-1)*3+1),3,1) + repmat([0;1;2],1,length(i));
            idx2 = repmat(((j-1)*3+1),3,1) + repmat([0;1;2],1,length(j));
            covMatrix = cov([data(:,idx1),data(:,idx2)]);
            switch lower(type)
                case 'cov'
                    tmp(j) = Escoufier_coefficient_fast(covMatrix);
                case 'cor'
                    tmp(j) = Escoufier_coefficient_fast(corrcov(covMatrix));
            end
        end
        RV(i,:) = tmp;
        parfor_progress;
    end
    closeParForProgress(path,ID);
    out = triu(RV,1)'+ triu(RV);
end