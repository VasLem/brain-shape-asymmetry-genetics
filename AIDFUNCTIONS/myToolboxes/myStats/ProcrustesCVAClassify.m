function out = ProcrustesCVAClassify(X,G,varargin)
         [n,nrV] = size(X);
         GL = unique(G);
         nrG = length(GL);
         results = zeros(n,2);
         % LOO validation
         index = (1:n);
         parfor i=1:n
             indexfor = setdiff(index,i);
             Xfor = X(indexfor,:);
             Gfor = G(indexfor);
             test = X(i,:);
             testG = G{i};
             testGnr = find(strcmp(testG,GL));
             [class,err] = classify(test,Xfor,Gfor);
             
             
         end
         





end