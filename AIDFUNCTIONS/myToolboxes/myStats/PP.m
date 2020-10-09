function [X,AvgA,AvgB] = PP(A,B)
         % get dimensions
            n = size(A,1);
            nA = size(A,2);
            nB = size(B,2);
            if ~n==size(B,1), error('A and B must have an equal number of rows (observations)'); end
         % Concatenate Data
            D = [B,A];
         % Center Data
            AvgD = mean(D,1);
            AvgA = AvgD(end-nA+1:end);AvgB = AvgD(1:nB);
            D = D - repmat(AvgD,n,1);
         % Computer singular Value decomposition
            if n < (nA + nB)
                [EigVec,Tcoeff,EigVal] = my_princomp(D,'econ');
            else
                [EigVec,Tcoeff,EigVal] = my_princomp(D);
            end
            clear Tcoeff;
         % Reduce the eigenvector system to the rows affecting A
            U = EigVec((nB+1:end),:);
            W = diag(EigVal);
            C = (W*U'/(U*W*U'));
         % Converting back to the B space
            out = repmat(AvgD,nA,1)' + EigVec*C;
            X = out(1:nB,:)';
end