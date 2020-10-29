function [out]=buildCRmatrix(data,Type)
% INPUT:
%
% data: nbParticipants x nbLMs, where nbLMS = 3 * nPoints on the face 
% Type: 1 to take the normalized CR matrices; 0 to take only the
% un-normalized CR matrix

% OUTPUT:
% out: nbLMs x nbLMs (different types of matrix inversion)

    nbpts = size(data,2)/3;
    CR = nan*zeros(nbpts);
    nrC = 3;

    disp('BUILDING CR MATRIX');
    [path,ID] = setupParForProgress(nbpts);
    
    MM = ones(nrC,nrC);
    MM(find(eye(nrC,nrC))) = 0;
    
    parfor i=1:nbpts
        tmp = zeros(1,nbpts);
        for j=i:nbpts
            idx1 = repmat(((i-1)*nrC+1),nrC,1) + repmat([0;1;2],1,length(i));
            idx2 = repmat(((j-1)*nrC+1),nrC,1) + repmat([0;1;2],1,length(j));
            S = cov([data(:,idx1),data(:,idx2)]);
            S = corrcov(S);
            
            S11 = S(1:nrC,1:nrC).*MM;
            S22 = S(nrC+1:end,nrC+1:end).*MM;
            S12 = S(1:nrC,nrC+1:end);
            S21 = S(nrC+1:end,1:nrC);          
            numCR = trace(S12*S21);
            denomCR = sqrt(trace(S11*S11)*trace(S22*S22));
            tmp(j) = sqrt(numCR/denomCR);           
%             s11_0 = S(1:nrC,1:nrC);
%             s22_0 = S(nrC+1:end,nrC+1:end);
%             s11_0(1:nrC+1:nrC*nrC) = 0; % set diag of s11 to zero; first value : step : final value
%             s22_0(1:nrC+1:nrC*nrC) = 0; % set diag of s22 to zero; first value: step : final value
%             s12 = S(1:nrC,nrC+1:end);
%             s21 = S(nrC+1:end,1:nrC);
%             % the covariance ratio
%             %CR(i,j) = sqrt( sum(diag(s12*s21))/sqrt(sum(diag(s11_0*s11_0)).*sum(diag(s22_0*s22_0))) );
%             tmp(j) = sqrt( sum(diag(s12*s21))/sqrt(sum(diag(s11_0*s11_0)).*sum(diag(s22_0*s22_0))) );
        end
        CR(i,:) = tmp;
        parfor_progress;
    end
    closeParForProgress(path,ID);
    
    out.CR = triu(CR,1)'+ triu(CR);
    
    if Type==1
%         outRowNorm = nan*zeros(nbpts);
%         for i=1:nbpts
%             outRowNorm(i,:) = out.CR(i,:)/max(out.CR(i,:)); % normalization 
%         end
        
        maxRow = max(out.CR,[],2);
        outRowNorm = out.CR./repmat(maxRow,1,size(out.CR,2));
         
        out.RN1 = triu(outRowNorm,1)'+triu(outRowNorm); % tril and triu matrices are different
        out.RN2 = tril(outRowNorm,-1)+tril(outRowNorm)';

        out.RN5 = bsxfun(@rdivide, out.CR, sqrt(sum(out.CR.^2, 2)));
        out.RN5 = triu(out.RN5,1)'+triu(out.RN5);

        out.RN7 = bsxfun(@rdivide, out.CR, sqrt(sum(out.CR.^2, 1)));
        out.RN7 = triu(out.RN7,1)'+triu(out.RN7);
    end
end