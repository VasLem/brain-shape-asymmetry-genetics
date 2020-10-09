function out = ShapeModelCVAClassify(X,G,AM,percvar,type,varargin)
         [n,~] = size(X);
         disp('Building Shape Model');
         SM = shapePCA;
         SM.RefScan = AM;
         getAverage(SM,X');
         getModel(SM,X');
         stripPercVar(SM,percvar);
         if strcmpi(type,'mahalanobis')
             disp('Normalizing');
             X = SM.Tcoeff./repmat(SM.EigVal',n,1);
         else
             X = SM.Tcoeff;
         end
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
             class = classify(test,Xfor,Gfor);
             resultfor = zeros(1,2);
             resultfor(2) = testGnr;
             resultfor(1) = double(strcmp(class,testG));
             results(i,:) = resultfor;
         end
         out.Correct = (sum(results(:,1))/n)*100;
         for i=1:1:nrG
            Group(i).Name = GL{i}; %#ok<*AGROW>
            index = find(results(:,2)==i);
            Group(i).Correct = (sum(results(index,1))/length(index))*100;
         end
         out.Group = Group;
         out.Results = results;
end