function Xnew = matrixInterShuffle(X)
         nX = size(X,1);
         nV = size(X,2);
         nL = size(X,3);
         Xnew = zeros(nL*nX,nV);
         ind = (1:nL*nX);
         for i=1:1:nL
             tmp = find(mod(ind,nL)==i-1);
             Xnew(tmp,:) = X(:,:,i); %#ok<FNDSB>
         end
end