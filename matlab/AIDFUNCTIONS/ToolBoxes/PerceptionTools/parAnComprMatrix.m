function out = parAnComprMatrix(CM,t)
         n = numel(CM);
         nr = size(CM,1);
         [pc,lc]=eigen(CM);
         F=pc*diag(sqrt(lc));
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
            indexfor = randsample(n,n,true);
            indexdiagfor = randsample(length(diagindex),length(diagindex),true);
            CMfor = CM;
            CMfor(index) = CM(indexfor);
            CMfor = (CMfor + CMfor')/2;
            tmp = diag(CM);
            CMfor(diagindex) = tmp(indexdiagfor);
            %CMfor(diagindex) = CM(diagindex);
            %CMfor(diagindex) = ones(1,length(diagindex));
            [~,lcfor]=eigen(CMfor);
            lcindex = 1:length(lcfor);
            dummyfor = dummy;
            dummyfor(lcindex) = lcfor;
            lcCount(:,i) = dummyfor>=lc;
         end
         toc;
         out.p = sum(lcCount,2)/t;
end