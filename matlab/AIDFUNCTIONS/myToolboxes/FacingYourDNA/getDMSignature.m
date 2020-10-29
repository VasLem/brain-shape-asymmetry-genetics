function [SingVal, SingVect]=getDMSignature(LM,k,ref)
% input:
% SymHead : 21450 x nbParticipant
% pRef: participant of reference
% p: new participant
% idx: index position for LMs
% k: nb of singular values

% output:
% SingVal: singular values sorted
% SingVect: singular vectors sorted
    if nargin<3, ref = [];end

    d=squareform(pdist(LM'));
    [u,s]=svds(d,k); % default 6 largest eigenvalues
    s = diag(s);

    if isempty(ref), SingVal = s; SingVect = u; return; end
    
    D = nan*zeros(k,k);
    for p=1:k
        for q=p:k
            D(p,q)=(dot(ref(:,p),u(:,q))/(norm(ref(:,p))*norm(u(:,q)))); % get angle between the reference vector and the new vector
        end
    end
    D = triu(D,1)'+ triu(D);
    idx = zeros(1,k);
    %signs = zeros(1,k);
    for j=1:k
        [~,idx(j)] = max(abs(D(j,:)));
        if sign(D(j,idx(j)))==-1, u(:,j) = -u(:,j);end
    end
    SingVal = s(idx);
    SingVect = u(:,idx);
end


% for j=1:k
%     if D(j,j)==-1
%         SingVect(:,j)=-u(:,j);
%     end
%     if abs(D(j,j))~=1
%         idx = find(abs(D(j,:))==max(abs(D(j,:)))); 
%         SingVect(:,j) = u(:,idx);
%         SingVal(j) = s(idx);
%     else
%         SingVect(:,j) = u(:,j);
%         SingVal(j) = s(j);
%     end
% end