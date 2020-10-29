function [SingVal, SingVect, D]=getSignature(SymHead,pRef,p,idx,k)
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

r=reshape(SymHead(:,p),3,7150);
r=r(:,idx);
d=squareform(pdist(r'));
[u,s,~]=svds(d,k); % default 6 largest eigenvalues
s = diag(s);

SingVect = nan*zeros(size(u));
SingVal = nan*zeros(size(s));


D = nan*zeros(k,k);
for p=1:k
    for q=p:k
        D(p,q)=round(dot(RefSingVect(:,p),u(:,q))/(norm(RefSingVect(:,p))*norm(u(:,q)))); % get angle between the reference vector and the new vector
    end
end

D = triu(D,1)'+ triu(D);

for j=1:k
    if D(j,j)==-1
        SingVect(:,j)=-u(:,j);
    end
    
    
    if abs(D(j,j))~=1
        idx = find(abs(D(j,:))==max(abs(D(j,:)))); 
        SingVect(:,j) = u(:,idx);
        SingVal(j) = s(idx);
    else
        SingVect(:,j) = u(:,j);
        SingVal(j) = s(j);
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