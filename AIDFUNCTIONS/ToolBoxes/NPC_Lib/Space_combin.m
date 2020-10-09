function [Space] = Space_combin(n1,n2,B)
%To generate the indexes of the permutation sample space.
%Space = Space_combin(n1,n2) generates the matrix NC-by-n1+n2 of all
%combinations of the n1+n2 observations with n1 elements. NC is equal to
%(n1+n2)!/n1![(n1+n2)-n1]!. In the last row the 1:(n1+n2) vector which
%refers to the observed sample.
%
%Space = Space_combin(n1,n2,B) to have only B samples of the permutation
%sample space.
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin==2
   B=0;
end

n=n1+n2;


if n>20&B==0
    B=iNPut('\n It can take a considerable amount of time, \n do you want to sample the space instead of \n comuputing the entire one? \n (ENTER=NO, scalar= Number of random permutations)')
    if isempty(B)==1;
        B=0;
    end
end

if B==0;
    Space=nchoosek(1:n,n1);
    %temp2=repmat(1:n,size(temp,1),1);
    
    for i=size(Space,1):-1:1
        Space(i,n1+(1:n2))=setdiff(1:n,Space(i,:));
    end
%     temp2=sort(temp2,2);
%     temp2=temp2(:,n1+1:n);
%     Space=[temp temp2];
Space=Space(end:-1:1,:);

else
    Space=zeros(B+1,n);
    for i=1:B
        Space(i,:)=randperm(n);
    end
    Space(B+1,:)=1:n;
end