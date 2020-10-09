%% Weighted Least Squares Superimposition using SVD

function T = superimpose_wLS_SVD(obj1, obj2, w)

% Function that superimposes the second object onto the first using 
% weighted least squares with SVD for the calculation of the rotation 
% matrix. This method expects that the first and the second object have 
% the same size.
% input:    objects obj1 and obj2 in homogeneous coordinates
%           w = vector or matrix with the weights
% output:   T = transformation matrix that superimposes the second object
%               onto the first

nbpts = size(obj1,2);
obj1 = obj1(1:3,:);     % no longer homogeneous coordinates
obj2 = obj2(1:3,:);

totalw = sum(w,2);
not0weight = find(sum(w,1));    % remove point having weight 0
if length(not0weight)~=nbpts
   obj1 = obj1(:,not0weight);
   obj2 = obj2(:,not0weight);
   w = w(:,not0weight);
   totalw = sum(w,2);
   nbpts = length(not0weight);
end
rw = repmat(w,3,1);

% Superimpose means by translating both objects to origin
centrum1 = sum(rw.*obj1,2)./repmat(totalw,3,1);
obj1 = obj1 - repmat(centrum1,1,nbpts);
centrum2 = sum(rw.*obj2,2)./repmat(totalw,3,1);
obj2 = obj2 - repmat(centrum2,1,nbpts);

% Scale both by making mean length 1
lengths = sqrt(sum((rw.*obj1).^2,1));
meanLength1 = sum(lengths,2)/totalw;
obj1 = obj1./meanLength1;
lengths = sqrt(sum((rw.*obj2).^2,1));
meanLength2 = sum(lengths,2)/totalw;
obj2 = obj2./meanLength2;
scalefactor = meanLength1/meanLength2;

% Rotate
[U,S,V] = svd(rw.*obj1*obj2');
H = V*sign(S)*U';
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