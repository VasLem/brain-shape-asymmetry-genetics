function [out]=buildCRmatrixOriginal(data,Type)
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

h=waitbar(0,'please wait');
for i=1:nbpts
    for j=i:nbpts
        idx1 = repmat(((i-1)*nrC+1),nrC,1) + repmat([0;1;2],1,length(i));
        idx2 = repmat(((j-1)*nrC+1),nrC,1) + repmat([0;1;2],1,length(j));
        S = cov([data(:,idx1),data(:,idx2)]);
        s11_0 = S(1:nrC,1:nrC);
        s22_0 = S(nrC+1:end,nrC+1:end);
        s11_0(1:nrC+1:nrC*nrC) = 0; % set diag of s11 to zero; first value : step : final value
        s22_0(1:nrC+1:nrC*nrC) = 0; % set diag of s22 to zero; first value: step : final value
        s12 = S(1:nrC,nrC+1:end);
        s21 = S(nrC+1:end,1:nrC);
        % the covariance ratio
        CR(i,j) = sqrt( sum(diag(s12*s21))/sqrt(sum(diag(s11_0*s11_0)).*sum(diag(s22_0*s22_0))) );
    end
    waitbar(i/nbpts,h);
end
close(h)

if Type==1
    out.CR = triu(CR,1)'+ triu(CR);
    outRowNorm = nan*zeros(nbpts);
    for i=1:nbpts
        outRowNorm(i,:) = out.CR(i,:)/max(out.CR(i,:)); % normalization 
    end
    
    out.RN1 = triu(outRowNorm,1)'+triu(outRowNorm); % tril and triu matrices are different
    out.RN2 = tril(outRowNorm,-1)+tril(outRowNorm)';
    
    out.RN5 = bsxfun(@rdivide, out.CR, sqrt(sum(out.CR.^2, 2)));
    out.RN5 = triu(out.RN5,1)'+triu(out.RN5);
    
    out.RN7 = bsxfun(@rdivide, out.CR, sqrt(sum(out.CR.^2, 1)));
    out.RN7 = triu(out.RN7,1)'+triu(out.RN7);
else
    out = triu(CR,1)'+ triu(CR);
end
end