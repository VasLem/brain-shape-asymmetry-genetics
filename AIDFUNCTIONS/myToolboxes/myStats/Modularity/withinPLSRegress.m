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