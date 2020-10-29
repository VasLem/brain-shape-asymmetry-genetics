function R2 = getLDStructure(GENO1,GENO2)
         pw = false;
         if nargin<2,GENO2 = GENO1; pw = true; end
         nG1 = size(GENO1,1);
         nG2 = size(GENO2,1);
         R2 = nan*zeros(nG1,nG2);
         %[path,ID] = setupParForProgress(nG1);tic;
         parfor i=1:nG1
             %i=607;
             site1 = double(GENO1(i,:));
             index1 = find(site1>=0);
             if isempty(index1), continue; end
             if length(unique(site1(index1)))==1, continue; end
             tmp = nan*zeros(1,nG2);
             for j=1:nG2
                if pw&&j<=i, continue; end% pairwise self computation not need to do all
                %j=862;
                site2 = double(GENO2(j,:));
                index2 = find(site2>=0);
                if isempty(index2), continue; end
                if length(unique(site2(index2)))==1, continue; end
                ind = intersect(index1,index2);
                if isempty(ind), continue; end
                out = mySNPLDEM(site1(ind),site2(ind));
                tmp(j) = out.r2;
             end
             R2(i,:) = tmp;
             %parfor_progress;
         end
         %closeParForProgress(path,ID);toc;
         if pw, R2(isnan(R2)) = 0;R2 = R2+R2'+eye(nG1,nG1); end
end