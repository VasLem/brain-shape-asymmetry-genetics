function [Space] = Space_sign(n,B)
%To Generate the +1 and -1  of the permutation sample space
%Space = Space_sign(n) generates the 2^n-by-n matrix of all permutations of
%n plus and minus signs. In the last row are the signs of the observed sample.
%
%Space = Space_sign(n,B) to have only B samples of the permutation
%sample space.
%
%This function is part of NPClib.
%NPC Lib Matlab library (companion software of F. Pesarin, L. Salmaso: Permutation Tests for Complex Problems: theory, applications and software, John Wiley and sons).
%developed by Livio Finos; consulting team: Rosa Arboretti, Francesco Bertoluzzo, Stefano Bonnini, Chiara Brombin, Livio Corain, Luigi Salmaso and Aldo Solari.
%E-mail for technical questions: livio@stat.unipd.it (Livio Finos).

if nargin==1
   B=0;
end

if B==0;
    Space=[-1;1];
    for i=2:n
        Space=[ Space -ones(2^(i-1),1); Space ones(2^(i-1),1)];
    end
else
    Space=[2.*binornd(1,.5,B,n)-1;ones(1,n)];
end