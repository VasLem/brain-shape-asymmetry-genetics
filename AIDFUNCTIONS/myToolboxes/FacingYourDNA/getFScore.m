function [S,SNPS,COVS,GBS] = getFScore(SNPM,SNPW,COVM,COVW,GBM,GBW)
        % Check input arguments
         if nargin < 2, SNPW = []; end % SNP weights
        % Initialize
         SNPS = 0;COVS = 0;GBS = 0;
         SNPN = 0;COVN = 0;GBN = 0;
        % SNP MATCHES 
         if ~isempty(SNPM)
             if isempty(SNPW)
                SNPN = size(SNPM,1);
                SNPW = ones(size(SNPM));
             else
                SNPN = nansum(SNPW);
                SNPW = repmat(SNPW,1,size(SNPM,2));
             end
             SNPS = nansum(SNPW.*-log(SNPM),1);
         end
         if nargin < 3, S = SNPS./SNPN; return; end
        % COV MATCHES 
         if ~isempty(COVM)
             if isempty(COVW)
                COVN = size(COVM,1); 
                COVW = ones(size(COVM));
             else
                COVN = nansum(COVW);
                COVW = repmat(COVW,1,size(COVM,2));
             end
             COVS = nansum(COVW.*-log(COVM),1);
         end
         if nargin < 5, S = (SNPS+COVS)./(SNPN+COVN); return; end
        % GB MATCHES
         if ~isempty(GBM)
             if isempty(GBW)
                GBN = size(GBM,1); 
                GBW = ones(size(GBM));
             else
                GBN = nansum(GBW); 
                GBW = repmat(GBW,1,size(GBM,2));
             end
             GBS = nansum(GBW.*-log(GBM),1);
         end
         S = (SNPS+COVS+GBS)./(SNPN+COVN+GBN);
end