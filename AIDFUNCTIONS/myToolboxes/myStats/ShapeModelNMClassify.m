function out = ShapeModelNMClassify(X1,X2,AM,percvar,type)
         %Warning Off
         [n1,~] = size(X1);
         [n2,~] = size(X2);
         n = n1+n2;
         disp('Building Shape Model');
         SM = shapePCA;
         SM.RefScan = AM;
         getAverage(SM,[X1' X2']);
         getModel(SM,[X1' X2']);
         stripPercVar(SM,percvar);
         X1 = SM.Tcoeff(1:n1,:);
         X2 = SM.Tcoeff(n1+1:end,:);
         results = zeros(n,3);
         resultsPP = zeros(n,3);
         resultsPar = zeros(n,3);
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
             D = getDirection(SM,AvgX1for,AvgX2for);
             if strcmpi(type,'mahalanobis')             
                [out1,outPP1,outPar1,angle1] = getDistance(SM,test,AvgX1for,'mahalanobis',D);
                [out2,outPP2,outPar2,angle2] = getDistance(SM,test,AvgX2for,'mahalanobis',D);
             else
                [out1,outPP1,outPar1,angle1] = getDistance(SM,test,AvgX1for,'euclidean',D);
                [out2,outPP2,outPar2,angle2] = getDistance(SM,test,AvgX2for,'euclidean',D);
             end
             % Original distance
             diff = out1-out2;
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
             % Perpendicular distance
             diff = outPP1-outPP2;
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
             resultsPP(i,:) = resultfor;
             % Parallell distance
             diff = outPar1-outPar2;
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
             resultsPar(i,:) = resultfor;
         end
         out.Correct = (sum(results(:,1))/n)*100;
         out.CorrectPP = (sum(resultsPP(:,1))/n)*100;
         out.CorrectPar = (sum(resultsPar(:,1))/n)*100;
         for i=1:1:2
            index = find(results(:,3)==i);
            Group(i).Correct = (sum(results(index,1))/length(index))*100;
            Group(i).Differences = results(index,2); %#ok<*AGROW>
            Group(i).CorrectPP = (sum(resultsPP(index,1))/length(index))*100;
            Group(i).DifferencesPP = resultsPP(index,2); %#ok<*AGROW>
            Group(i).CorrectPar = (sum(resultsPar(index,1))/length(index))*100;
            Group(i).DifferencesPar = resultsPar(index,2); %#ok<*AGROW>
         end
         out.Group = Group;
         out.Results = results;
         out.ResultsPP = resultsPP;
         out.ResultsPar = resultsPar;
end