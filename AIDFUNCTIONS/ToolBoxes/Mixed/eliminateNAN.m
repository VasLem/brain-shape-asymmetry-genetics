function [A,B,index] = eliminateNAN(A,B)
         index = (1:size(A,1));
         [i,~] = find(isnan(A));
         i = unique(i);
         index = setdiff(index,i);
         A = A(index,:);
         B = B(index,:);
end