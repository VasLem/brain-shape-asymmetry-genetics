%% Least Squares Superimposition using SVD

function T = superimpose_LS_SVD(obj1, obj2)

% Function that superimposes the second object onto the first using least
% squares with SVD for the calculation of the rotation matrix. This method
% expects that the first and the second object have the same size.
% input:    object in homogeneous coordinates
% output:   T = transformation matrix that superimposes the second object
%               onto the first

dim = size(obj1,1)-1;
nbpts = size(obj1,2);
obj1 = obj1(1:3,:);     % no longer homogeneous coordinates
obj2 = obj2(1:3,:);

% Superimpose means by translating both objects to origin
centrum1 = sum(obj1,2)./nbpts;
obj1b = obj1 - repmat(centrum1,1,nbpts);
centrum2 = sum(obj2,2)./nbpts;
obj2b = obj2 - repmat(centrum2,1,nbpts);

% Scale both by making mean length 1
lengths = sqrt(sum(obj1b.^2,1));
meanLength1 = sum(lengths,2)/nbpts;
obj1b = obj1b./meanLength1;
lengths = sqrt(sum(obj2b.^2,1));
meanLength2 = sum(lengths,2)/nbpts;
obj2b = obj2b./meanLength2;
scalefactor = meanLength1/meanLength2;

% Rotate
[U,S,V] = svd(obj1b*obj2b');
S = diag(S);
S = diag(sign(S));
H = V*S*U';
H = H';

% Put together transformation matrix
% First translate second object back to the origin
T1 = eye(4);
T1(1:3,4) = -1.*centrum2;
% Next scale the object
T2 = eye(4);
T2(1:3,1:3) = scalefactor*eye(3);
% Then rotate
T3 = eye(4);
T3(1:3,1:3) = H;
% Finally translate to the centre of the first object
T4 = eye(4);
T4(1:3,4) = centrum1;
% The final transformation matrix is the combination of all of this
T = T4*T3*T2*T1;
