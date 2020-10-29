function out = NonParametricManovaAsymmetry(X1,X2,t)
% I = Individuals
% DA = Directional Asymmetry
% FA = Fluctuating Asymmetry
% E = Error
         if nargin < 3, t = 0; end
         X1S = matrixInterShuffle(X1);
         X2S = matrixInterShuffle(X2);        
         %Following the nomenclature from M.J. Anderson, A new method for
         %non-parametric multivariate analysis of variance
         a = 2;
         b = size(X1,1);
         n = size(X1,3);
         % Building the Distance Matrix
         D = squareform(pdist([X1S;X2S],'euclidean'));
         % extraction Mean Sums of Squares
         SS = extractSS(D,a,b,n);
         % Generating Output
         out.MSI = SS.MSB;% Mean Sums of squares I
         out.FI = SS.MSB/SS.MSAB;% F-ratio I against FA
         out.MSD = SS.MSA;% Mean Sums of squares DA
         out.FD = SS.MSA/SS.MSAB;% F-ratio DA against FA
         out.MSF = SS.MSAB;% Mean Sums of squares FA
         out.FF = SS.MSAB/SS.MSR;% F-ratio of FA against E
         out.MSE = SS.MSR;      
         % Disp, permuting test
         if t==0, return; end
         disp('Permuting');
         tic; 
         FICount = false(1,t);
         FDCount = false(1,t);
         FFCount = false(1,t);
         avgX1 = mean(X1S);
         avgX2 = mean(X2S);
         avgR = (X1S+X2S)/2;
         avgG = mean([X1S;X2S]);
         parfor i=1:t
             % Colom-wise shuffeling of cells for Directional effect
                    Set1 = X1S;
                    Set2 = X2S;
                    r = randi(2,n*b,1);
                    index = find(r==2);
                    Set1(index,:) = X2S(index,:);
                    Set2(index,:) = X1S(index,:,:);
                    Dperm = squareform(pdist([Set1;Set2],'euclidean'));
                    SSperm = extractSS(Dperm,a,b,n);
                    FDperm = SSperm.MSA/SSperm.MSAB;% F-ratio DA against FA
                    FDCount(i) = FDperm >= out.FD;
                % Row-wise shuffeling for Individual effect
                    index = randperm(b);
                    Set1 = X1(index,:,:); %#ok<*PFBNS>
                    Set1 = matrixInterShuffle(Set1);
                    Dperm = squareform(pdist([Set1;X2S],'euclidean'));
                    SSperm = extractSS(Dperm,a,b,n);
                    FIperm = SSperm.MSB/SSperm.MSAB;% F-ratio I against FA
                    FICount(i) = FIperm >= out.FI;
                % Residual shuffeling for Interaction effect
                    Set1 = X1S-repmat(avgX1,n*b,1)-avgR+repmat(avgG,n*b,1);
                    Set2 = X2S-repmat(avgX2,n*b,1)-avgR+repmat(avgG,n*b,1);
                    Set = [Set1;Set2];
                    index = randperm(n*b*2);
                    Set = Set(index,:);
                    Set1 = Set(1:n*b,:);
                    Set2 = Set(n*b+1:end,:);
                    Dperm = squareform(pdist([Set1;Set2],'euclidean'));
                    SSperm = extractSS(Dperm,a,b,n);
                    FFperm = SSperm.MSAB/SSperm.MSR;% F-ratio of FA against E
                    FFCount(i) = FFperm>=out.FF;
         end
         toc;
         out.FIPperm = sum(FICount)/t;
         out.FDPperm = sum(FDCount)/t;
         out.FFPperm = sum(FFCount)/t;
end
function Xnew = matrixInterShuffle(X)
         nX = size(X,1);
         nV = size(X,2);
         nL = size(X,3);
         Xnew = zeros(nL*nX,nV);
         ind = (1:nL*nX);
         for i=1:1:nL
             tmp = find(mod(ind,nL)==i-1);
             Xnew(tmp,:) = X(:,:,i); %#ok<FNDSB>
         end
end
function out = extractSS(D,a,b,n)
         N = a*b*n;
         T = 0;A = 0;B = 0;R = 0;
         for i=1:1:(N-1)
             for j=(i+1):1:N
                 d = D(i,j)^2;
                 T = T+d;
                 if inA(i,N,a)==inA(j,N,a)
                    A = A+d;
                 end
                 if inB(i,N,a,b)==inB(j,N,a,b)
                    B = B+d;
                    if inA(i,N,a)==inA(j,N,a)
                       R = R+d;
                    end
                 end
             end
         end
         out.SST = T/N;
         out.SSAW = A/(n*b);
         out.SSA = out.SST-out.SSAW;
         out.SSBW = B/(n*a);
         out.SSB = out.SST-out.SSBW;
         out.SSR = R/n;
         out.SSAB = out.SST-out.SSA-out.SSB-out.SSR;
         out.MSA = out.SSA/(a-1);
         out.MSB = out.SSB/(b-1);
         out.MSAB = out.SSAB/((a-1)*(b-1));
         out.MSR = out.SSR/(N-a*b);
end
function out = inA(i,N,a)
         nA = N/a;
         out = ceil(i/nA);
end
function out = inB(i,N,a,b)
         nA = N/a;
         nB = (N/a)/b;
         i = i-((inA(i,N,a)-1)*nA);
         out = ceil(i/nB);
end

