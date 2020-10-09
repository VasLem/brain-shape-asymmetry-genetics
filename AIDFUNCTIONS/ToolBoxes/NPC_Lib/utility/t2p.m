function[P T]=t2p(T,Pobs,tail)
%P=t2p(T,Pobs,tail)
%Compute the conjoint p-values permutation space.
%
% [P]=t2p(T) return a P matrix (the conjoint p-values permutation space) 
% using the matrix of statistics T.
% [P]=t2p(T,1) return a P vector of observed p-values 
% using the matrix of statistics T.
% [P]=t2p(T,0) is the same as  [P]=t2p(T).
%
%
%References: Pesarin, F. (2001) Multivariate Permutation Test 
%             with Application in Biostatistics. Wiley, New York.
%
%Livio Finos
%e-mail: livio@stat.unipd.it

if isempty(T)
    P=T;
else
sizorig=size(T);
T=reshape(T,sizorig(1),prod(sizorig(2:end)));

[B m]=size(T);
if nargin==3
    if prod(size(tail))==1
        tail=repmat(tail,[1 m]);
    else
        if not(prod(size(tail))==prod(sizorig(2:end))) 
            fprintf('\nWarning: tail matrix does not match the size of T; the first value tail(1) will be applied to all other variables'); 
            tail=repmat(tail(1),[1 m]);
        end
    end

        T(:,find(tail==0))=abs(T(:,find(tail==0)));
        T(:,find(tail==-1))=-T(:,find(tail==-1));
end

%T=single(T);

if nargin==1
    Pobs=0;
end

sizorigP=sizorig;
if Pobs==1
    for j=1:m
        P(1,j)=sum(T(end,j)<=T(:,j));%+sum(T(end,j)==T(:,j))./2;
    end
    sizorigP(1)=1;
else
    P=maxrank(-T);
end


P=P./B;
P=reshape(P,sizorigP);
T=reshape(T,sizorig);
end


