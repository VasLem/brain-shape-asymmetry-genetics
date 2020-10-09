function [R,tS,S] = getRank(TM,FM,wM)
         if nargin < 9, GBW = []; end
         if nargin < 8, GBM = []; end
         if nargin < 7, tGBM = []; end
         if nargin < 6, COVW = []; end
         if nargin < 5, COVM = []; end
         if nargin < 4, tCOVM = []; end
         if nargin < 3, SNPW = []; end
         tS = getFScore(tSNPM,SNPW,tCOVM,COVW,tGBM,GBW);
         S = getFScore(SNPM,SNPW,COVM,COVW,GBM,GBW);
         R = sum(S<tS)+1;
%          SNPM = rGMatches(:,:,t);
%          SNPW = freqSNPs(:,t);
%          COVM = rSMatches(1,:,t); 
%          COVW = 1;
%          GBM = [];
%          GBW = [];      
%          N = nansum(SNPW) + 1 + sum(GBW);
%          testScore = (nansum(SNPW.*-log(testSNP)) + ...
%                      -log(testS) + ...
%                      nansum(GBW.*-log(testGB)))./N;
%          Score = (nansum(repmat(SNPW,1,size(SNP,2)).*-log(SNP)) + ...
%                  -log(S) + ...
%                  nansum(repmat(GBW,1,size(GB,2)).*-log(GB)))./N;             
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
%          rank = sum(Score<testScore)+1;
end