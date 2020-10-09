function out = ProcustesNMClassify(X1,X2,W)
         %Warning Off
         if nargin < 3, W = []; end
         [n1,~] = size(X1);
         [n2,~] = size(X2);
         n = n1+n2;
         results = zeros(n,3);
         % LOO validation
         parfor i=1:n
             switch i<=n1
                 case 1
                     indexfor = setdiff((1:n1),i);
                     X1for = X1(indexfor,:); %#ok<*PFBNS>
                     test = X1(i,:);
                     testG = 1;
                     X2for = X2;         
                 case 0
                     indexfor = setdiff((1:n2),i-n1);
                     X2for = X2(indexfor,:);
                     test = X2(i-n1,:);
                     testG = 2;
                     X1for = X1;
             end
             % constructing averages
             AvgX1for = mean(X1for);
             AvgX2for = mean(X2for);
             out1 = getProcrustesDistance(test,AvgX1for,W);
             out2 = getProcrustesDistance(test,AvgX2for,W);
             if isempty(W)
                diff = out1.RMSE-out2.RMSE;
             else
                diff = out1.WRMSE-out2.WRMSE;
             end
             resultfor = zeros(1,3);
             resultfor(3) = testG;
             resultfor(2) = diff;
             switch testG
                 case 1
                    if diff<=0
                       resultfor(1) = 1;
                    else
                       resultfor(1) = 0;
                    end
                 case 2
                    if diff<=0
                       resultfor(1) = 0;
                    else
                       resultfor(1) = 1;
                    end
             end
             results(i,:) = resultfor;
         end
         out.Correct = (sum(results(:,1))/n)*100;
         for i=1:1:2
            index = find(results(:,3)==i);
            Group(i).Correct = (sum(results(index,1))/length(index))*100;
            Group(i).Differences = results(index,2); %#ok<*AGROW>
         end
         out.Group = Group;
         out.Results = results;
end