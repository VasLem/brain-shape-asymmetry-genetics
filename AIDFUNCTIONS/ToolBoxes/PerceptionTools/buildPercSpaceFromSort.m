function out = buildPercSpaceFromSort(sort,t)
         [no,nj] = size(sort);
         % converting sorting task into Indicator, co-occurence and distance matrices
         [L,R,D] = Sort2LRD(sort);
         % Convert Distance matrices into standard cross product matrices
         [rS,S] = Dist2StandCrossProd(D);
         % Compute the inter assosor RV coefficients and extract assessor weights
         [RV] = getRvCoeff(S);
         [W] = Rv2Weights(RV,true);
         % Compute comprimise matrix
         C = getCompromise(S,W);
         % get eigen decomposition and F scores
         [pc,lc]=eigen((C+C')/2);
         F=pc*diag(sqrt(lc));       
         out.sort = sort;
         out.L = L;
         out.R = R;
         out.D = D;
         out.S = S;
         out.RV = RV;
         out.W = W;
         out.C = C;
         out.pc = pc;
         out.lc = lc;
         out.F = F;
         if t < 1, return; end
         lcCount = false(length(lc),t);
         nrlc = length(lc);
         tic;
         parfor i=1:t
             succes = 0;
             while ~succes
                 % permuting the sorting matrix
                 sortfor = sort;
                 for j=1:nj
                     ind = randperm(no);
                     sortfor(:,j) = sort(ind,j);
                 end
                 [~,~,Dfor] = Sort2LRD(sortfor);
                 [~,Sfor] = Dist2StandCrossProd(Dfor);
                 RVfor = getRvCoeff(Sfor);
                 Wfor = Rv2Weights(RVfor,false);
                 Cfor = getCompromise(Sfor,Wfor);
                 [~,lcfor]=eigen((Cfor+Cfor')/2);
                 if ~length(lcfor)==nrlc
                    dummy = zeros(1,nrlc);
                    dummy(1:length(lcfor)) = lcfor;
                    lcfor = dummy;
                 end
                 try
                    lcCount(:,i) = lcfor>=lc;
                    succes = 1;
                 catch % the random generated input gave crap, (this happens sometimes)
                     succes = 0;
                 end
             end
         end
         toc;
         out.p = sum(lcCount,2)/t;
end