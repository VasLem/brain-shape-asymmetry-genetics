function out = parAnComprMatrixvRandom(CM,t)
         n = numel(CM);
         nr = size(CM,1);
         
%          tmp = eye(nr,nr);
%          diagindex = find(tmp);
%          CM(diagindex) = ones(1,length(diagindex));
         
         [pc,lc]=eigen(CM);
         F=pc*diag(sqrt(lc));
         mu = mean(CM);
         sigma = std(CM);
         
         tmp = eye(nr,nr);
         diagindex = find(tmp);
         mudiag = mean(CM(diagindex));
         sigmadiag = std(CM(diagindex));      

         out.pc = pc;
         out.lc = lc;
         out.F = F;
         if t<1, return; end
         lcCount = false(length(lc),t);
         index = (1:n);
         tmp = eye(nr,nr);
         diagindex = find(tmp);
         dummy = zeros(size(lc));
         tic;
         parfor i=1:t
%             indexfor = randsample(n,n,true);
%             indexdiagfor = randsample(length(diagindex),length(diagindex),true);
%             CMfor = CM;
%             CMfor(index) = CM(indexfor);
%             CMfor = (CMfor + CMfor')/2;
%             tmp = diag(CM);
%             CMfor(diagindex) = tmp(indexdiagfor);
            %CMfor(diagindex) = CM(diagindex);
            %CMfor(diagindex) = ones(1,length(diagindex));
            CMfor = repmat(mu,nr,1) + repmat(sigma,nr,1).*randn(nr,nr);
            CMfor = (CMfor + CMfor')/2;
            CMfor(diagindex) = repmat(mudiag,1,length(diagindex)) + sigmadiag.*randn(1,length(diagindex));
%             CMfor(diagindex) = ones(1,length(diagindex));
            [~,lcfor]=eigen(CMfor);
            lcindex = 1:length(lcfor);
            dummyfor = dummy;
            dummyfor(lcindex) = lcfor;
            lcCount(:,i) = dummyfor>=lc;
         end
         toc;
         out.p = sum(lcCount,2)/t;
end