function out = extractModules(X,nm)
         nV = size(X,2)/3;
         nO = size(X,1);
         avg = LMObj('Vertices',reshape(mean(X),3,nV));
         if nargin < 2,
            nm = floor(nO/3);
         end
         index = zeros(1,nm);
         varmap = zeros(nm,nV);
         XCopy = X;
         for i=1:1:nm
             % all PLS regressions
                pctvar = withinPLSRegress(XCopy);
             % find the most explaining point
                [val, index(i)] = max(pctvar);
             % perform the regression of that point on the rest
                A = [XCopy(:,(index(i)-1)*3+1) XCopy(:,(index(i)-1)*3+2) XCopy(:,(index(i)-1)*3+3)];
                [~,~,~,~,~,~,~,statR] = plsregress(A,XCopy,3);
             % work with the residual information
                XCopy = statR.Yresiduals;
                varmap(i,:) = pctvar;
         end
         out.index = index;
         out.varmap = varmap;
         % perform a partial regressions
         % budling points together
         A = [];% TO DO optimise programming
         for i=1:1:nm
             A = [A X(:,(index(i)-1)*3+1) X(:,(index(i)-1)*3+2) X(:,(index(i)-1)*3+3)];
         end
         Tind = (1:size(A,2));
         LocalR2 = zeros(nm,nV);
         for i=1:1:nm
                 ind = [(i-1)*3+1 (i-1)*3+2 (i-1)*3+3];
                 TRA = setdiff(Tind,ind);
             % Take Out effect of other centers on face points
                 [~,~,~,~,~,~,~,statsR] = plsregress(A(:,TRA),X);
                 E = statsR.Yresiduals;
             % Take Out effect of other centers on this center
                [~,~,~,~,~,pctvar,~,statsR] = plsregress(A(:,TRA),A(:,ind));
             % partial regression
                 Afor = statsR.Yresiduals;
                 [~,~,~,~,M,pctvar,~,stats] = plsregress(Afor,E,3);
                 [~,LocalR2(i,:)] = getLocalStatistic(M,stats.Yresiduals,E,'shape');
         end
         % Classification
         [maxval,maxind] = max(LocalR2);
         viewer(avg);
         avg.Value = maxind;
         avg.ColorMode = 'Indexed';
         out.LocalR2 = LocalR2;
         out.avg = avg;
         out.maxval = maxval;
         out.maxind = maxind;
         % perform a multiple regression
         MR = multiplePlsRegress(A,X,0,'shape');
         out.MR = MR;
end

function [LocalE,LocalS] = getLocalStatistic(M,R,B,type)
         [n,nB] = size(B); 
         P = B-R;
         A = repmat(mean(B),n,1);
         SST = sum((B-A).^2);
         SSR = sum((P-A).^2);
         switch lower(type)
             case 'value'
                 LocalE = M(2,:);
                 LocalS = SSR./SST;
              case 'shape'
                 LocalE = zeros(4,nB/3);
                 LocalE(1:3,:) = reshape(M(2,:),3,nB/3);
                 LocalE(4,:) = sqrt(sum(LocalE(1:3,:).^2));
                 SSR = reshape(SSR,3,nB/3);
                 SSR = sum(SSR);
                 SST = reshape(SST,3,nB/3);
                 SST = sum(SST);
                 LocalS = SSR./SST;
             otherwise
                 error('unknown type');
         end        

end

function [pctvar] = withinPLSRegress(X)
         nrV = size(X,2)/3;
         pctvar = zeros(1,nrV);
         tic;
         parfor i=1:nrV
                A = [X(:,(i-1)*3+1) X(:,(i-1)*3+2) X(:,(i-1)*3+1)];
                [~,~,~,~,~,forpctvar] = plsregress(A,X,3);
                pctvar(i) = sum(forpctvar(2,:));           
         end
         toc;
end
