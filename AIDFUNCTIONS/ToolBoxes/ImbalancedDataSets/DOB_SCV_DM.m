function F = DOB_SCV_DM(K,D,G)
    %if isscalar(D), D = randn(D,D);end
    n = size(D,1);
    if nargin < 3, G = ones(n,1); end% simply distribute similar data features accross
    if isempty(G), G = ones(n,1); end% one big group
    % dealing with unknown group memberships
    indnan = find(isnan(G)); %#ok<*EFIND>
    if isempty(indnan)
        F = nan*zeros(1,n);
        Fid = (2:K);
        lC = unique(G);
        nC = length(lC);
        Gfor = G;
        for i=1:nC
           indC = find(Gfor==lC(i));
           nrC = length(indC);
           while nrC>0
              e = randsample(1:nrC,1);
              F(indC(e)) = 1;% assign to first fold
              Gfor(indC(e)) = nan;% reducing the group members to deplete them and stop the while
              if nrC==1,break;end
              re = setdiff(1:nrC,e);
              distances = D(indC(e),indC(re));
              %distances = sqrt(sum((repmat(D(indC(e),:),length(re),1)-D(indC(re),:)).^2,2));
              [~,index] = sort(distances,'ascend');
              if length(index)<K-1
                 F(indC(re(index))) = Fid(1:length(index));
                 Gfor(indC(re(index))) = nan;
              else
                 % assign to folds
                 F(indC(re(index(1:K-1)))) = Fid;
                 Gfor(indC(re(index(1:K-1)))) = nan;
              end
              indC = find(Gfor==lC(i));
              nrC = length(indC);
           end    
        end
    else
       indrest = setdiff((1:n),indnan);
       Fnan = DOB_SCV_DM(K,D(indnan,indnan));% distribute nan values according to similarity in D only
       Frest = DOB_SCV_DM(K,D(indrest,indrest),G(indrest));% distribute as originally planned
       F(indnan) = Fnan;
       F(indrest) = Frest;
    end
end