function [Space] = Space_perm(N,B,connected)
%To generate the indexes of the permutation sample space.
%Space = Space_perm(N) generates the N!-by-sum(N) matrix indexes of all
%permutations of the N observations. If the observations are stratified
%N must be a vector with m elements where m is the number of strata. In the
%last row is 1:N vector which refers to the observed sample.
%
%Space = Space_perm(N,B) to have only B samples of the permutation sample
%space. If B = 0, it peeforms all permutations.
%
%Space = Space_perm(N,B,connected) if it is assumed the connection respects the
%observation ordering (i.e. first, second.. observations of each stratum)
%set connected different from 0. In this case the elements of N have to be
%equal. connected=0 by default.
% 
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).


if nargin==1
   B=0;
   connected=0;
elseif nargin==2
   connected=0;
end


if (max(N)>10)&(B==0)
    B=input('\n It can take a considerable amount of time, \n do you want to sample the space instead of \n computing the entire one? \n (ENTER=NO, scalar= Number of random permutations)')
    if isempty(B)==1;
        B=0;
    end
end

SPACE=[];
for j=1:length(N)
    n=N(j);
    if (j==1) | (connected==0)
        if B==0;
            Space=perms(1:n);
            Space=Space(end:-1:1,end:-1:1);
        else
            Space=zeros(B+1,n);
            for i=1:B
                Space(i,:)=randperm(n);
            end
            Space(B+1,:)=1:n;
        end
    end
    SPACE=[SPACE Space+size(SPACE,2)];
end
Space=SPACE;