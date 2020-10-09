% 2D Least Squares

function T = LS_2D(A,B)

% Function that superimposes the second object onto the first using least
% squares with SVD for the calculation of the rotation matrix. This method
% expects that the first and the second object have the same size.
% input:    object in homogeneous coordinates
% output:   T = transformation matrix that superimposes the second object
%               onto the first

dim = 2;
nbpts = size(A,2);
A = A(1:dim,:);     % no longer homogeneous coordinates
B = B(1:dim,:);

% Center the coordinqtes of one object on the other
centrumB = sum(B,2)./nbpts;
B = B - repmat(centrumB,1,nbpts);
centrumA = sum(A,2)./nbpts;
A = A - repmat(centrumA,1,nbpts);

% Scale both objects
lengths = sqrt(sum(A.^2,1));
meanLengthA = sum(lengths,2)/nbpts;
A = A./meanLengthA;
lengths = sqrt(sum(B.^2,1));
meanLengthB = sum(lengths,2)/nbpts;
B = B./meanLengthB;
scalefactor = meanLengthA/meanLengthB;

% Rotate
[U,S,V] = svd(A*B');
S = diag(S);
S = diag(sign(S));
H = V*S*U';
H = H';

% Make transformation matrix
% First translate second object back to the origin
T1 = eye(dim+1);
T1(1:dim,dim+1) = -centrumB;
% Next scale the object
T2 = eye(dim+1);
T2(1:dim,1:dim) = scalefactor*eye(2);
% Then rotate
T3 = eye(dim+1);
T3(1:dim,1:dim) = H;
% Finally translate to the centre of the first object
T4 = eye(dim+1);
T4(1:dim,dim+1) = centrumA;
% The final transformation matrix is the combination of all of this
T = T4*T3*T2*T1;