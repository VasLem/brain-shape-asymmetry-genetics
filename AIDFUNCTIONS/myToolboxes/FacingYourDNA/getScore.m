function [S,SNPS,COVS,GBS] = getScore(SNPM,SNPW,COVM,COVW,GBM,GBW)
        % Check input arguments
         if nargin < 2, SNPW = []; end % SNP weights
        % Initialize
         S = 0;SNPS = 0;COVS = 0;GBS = 0;
         N = 0;SNPN = 0;COVN = 0;GBN = 0;
        % SNP MATCHES 
         if ~isempty(SNPM)
             if isempty(SNPW)
                SNPW = ones(size(SNP));
                SNPN = size(SNP,1);
             else
                SNPW = repmat(SNPW,1,size(SNP,2));
                SNPN = nansum(SNPW);
             end
             SNPS = nansum(SNPW.*-log(SNPM));
         end
         if nargin < 3, S = SNPS./SNPN; return; end
        % COV MATCHES 
         if ~isempty(COVM)
             if isempty(SNPW)
                SNPW = ones(size(SNP));
                SNPN = size(SNP,1);
             else
                SNPW = repmat(SNPW,1,size(SNP,2));
                SNPN = nansum(SNPW);
             end
             SNPS = nansum(SNPW.*-log(SNPM));
         end
        % GB MATCHES
         
            
         
         
         
         
         
         if nargin < 6, GBW = []; end % GB weights
         if nargin < 5, GBM = []; end % GB matches
         if nargin < 4, COVW = []; end % COV weights
         if nargin < 3, COVM = []; end % COV matches
         
         
         
         
         N = nansum(SNPW) + 1 + sum(GBW);
         testScore = (nansum(SNPW.*-log(testSNP)) + ...
                     -log(testS) + ...
                     nansum(GBW.*-log(testGB)))./N;
         Score = (nansum(repmat(SNPW,1,size(SNP,2)).*-log(SNP)) + ...
                 -log(S) + ...
                 nansum(repmat(GBW,1,size(GB,2)).*-log(GB)))./N;             
%          testScore = 0;Score = 0;
%          nr = 0;
%          if ~isempty(testSNP)
%             if isempty(SNPW), SNPW = ones(size(testSNP)); end
%             %N = nansum(SNPW);
%             N = 1;
%             testScore = testScore + nansum(SNPW.*-log(testSNP))./N;
%             Score = Score + nansum(repmat(SNPW,1,size(SNP,2)).*-log(SNP))./N;
%             nr = nr +1;
%          end
%          if ~isempty(testS)
%             testScore = testScore + -log(testS);
%             Score = Score + -log(S);
%             nr = nr +1;
%          end
%          if ~isempty(testGB)
%             if isempty(GBW), GBW = ones(size(testGB)); end
%             %N = sum(GBW);
%             N = 1;
%             testScore = testScore + nansum(GBW.*-log(testGB))./N;
%             Score = Score + nansum(repmat(GBW,1,size(GB,2)).*-log(GB))./N;
%             nr = nr+1;
%          end
%          %Score = Score/nr;
%          %testScore = testScore/nr;
         rank = sum(Score<testScore)+1;
end