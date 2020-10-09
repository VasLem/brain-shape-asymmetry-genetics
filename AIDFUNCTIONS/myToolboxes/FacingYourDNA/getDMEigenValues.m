function [SingVal, SingVect]=getDMEigenValues(SymHead,pRef,p,idx,k)
% input:
% SymHead : 21450 x nbParticipant
% pRef: participant of reference
% p: new participant
% idx: index position for LMs
% k: nb of singular values

% output:
% SingVal: singular values sorted
% SingVect: singular vectors sorted


[~,RefSingVect]=getRefSignature(SymHead,pRef,idx,k);
nV = size(SymHead,1)/3;
r=reshape(SymHead(:,p),3,nV);
r=r(:,idx);
d=squareform(pdist(r'));
[u,s]=svds(d,k); % default 6 largest eigenvalues
s = diag(s);

SingVect = nan*zeros(size(u));
SingVal = nan*zeros(size(s));

for j=1:k
    disp(num2str(j))
    D=round(dot(u(:,j),RefSingVect(:,j))/(norm(RefSingVect(:,j))*norm(u(:,j)))); % get angle between the reference vector and the new vector
    
    if D==-1
        SingVect(:,j)=-u(:,j);
    end
    if abs(D)~=1 % singular vectors mismatch
               
        SingVect(:,j)=u(:,j+1);
        SingVect(:,j+1) = u(:,j);
        SingVal(p)=s(p+1);
        SingVal(p+1)=s(p);
    else % singular vectors match
        SingVect(:,j)=u(:,j);
        SingVal(j)=s(j);
    end
end

end



function [RefSingVal, RefSingVect]=getRefSignature(SymHead,pRef,idx,k)
% input:
% SymHead : 21450 x nbParticipant
% pRef: participant of reference
% idx: index position for LMs
% k: nb of singular values

% output:
% RefSingVal: reference singular values
% RefSingVect: reference singular vectors
r=reshape(SymHead(:,pRef),3,7150);
r=r(:,idx);
d=squareform(pdist(r'));
[RefSingVect,s,~]=svds(d,k);
RefSingVal = diag(s);
end