function X = tls(A,B)

% solves the TLS problem with possibly multiple RHS
% A X \approx B
% with A = m x n, B = m x d

% D. Sima, KULeuven, 2009

[m,n] = size(A);
[mm,d] = size(B);
if (m~=mm) 
    error('A and B should have the same number of rows');
end

% use QR before SVD for efficiency
R = qr([A B]); R = triu(R);
[U,Sigma,V] = svd(R);

% pad with zero in case [A B] is "fat"
Sigma = [diag(Sigma); zeros(n+1-m,1)]; 

% avoid non-existence of solution in nongeneric cases
k = n;
while (k>0) && ((Sigma(k)-Sigma(k+1)<eps) || (rank(V(n+1:n+d,k+1:end))~=d))
    k = k-1;
end

X = -V(1:n,k+1:end)/V(n+1:n+d,k+1:end);
