%% Weighted Least Squares Superimposition using iterative calculation

function T = superimpose_wLS_moment(obj1, obj2, w)

% Function that superimposes the second object onto the first using 
% weighted least squares with iterative calculation of the moment for the 
% calculation of the rotation matrix. This method expects that the first 
% and the second object have the same size.
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

% Rotate    %%aangepaste code van Peter%%
% Initialization
mintheta = 0.0001;
theta = pi/4;
H = eye(3);
% iterate rotation
totalerror = calcep(obj1,obj2,rw);
while theta > mintheta
    M = calcmoment(obj1,obj2,rw);
    if (M==0), break; end
    M = M/norm(M);
    rotation = rodrigues(M*theta);
    newObj2 = rotation*obj2;
    error = calcep(obj1,newObj2,rw);
    if (error<totalerror)
        obj2 = newObj2;
        H = rotation*H;
        totalerror = error;
    else
        theta = theta/2;
    end 
end

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

end

function ep = calcep(q,p,w)    %%aangepaste code van Peter%%
    ep = sum(sum(w.*((p-q).^2)));
end

function M = calcmoment(q,p,w)    %%aangepaste code van Peter%%
    M = sum(w.*cross(p,(q-p)),2);
end