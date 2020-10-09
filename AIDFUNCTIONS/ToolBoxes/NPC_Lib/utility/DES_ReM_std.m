function[DES]=DES_ReM_std(C,TYPE)  
% DES_ReM_std(C,TYPE) create a desing matrix for Repeated measures analysis
% C: number of samples
% TYPE:
%   Seq (default): sequential comparisons (1vs2, 2vs3, etc.)
%   All: All comparisons
%   Bal: Base line vs others (Base line is the first one)
%   Trd: Trend (1vs2...C, 12vs3...C, 123vs4...C)
%
%Livio Finos
%e-mail: livio@stat.unipd.it

if nargin==1
    TYPE='seq';
else
	TYPE=lower(TYPE);
end

switch TYPE
    case {'seq'}
        DES=-eye(C);DES=DES(:,1:end-1);
        piu=find(DES==-1)+1;
        DES(piu)=1;
    case {'all'}
        DES=zeros(C,0);
        for i=1:C-1
        D=-eye(C);D=D(:,1:end-i);
        piu=find(D==-1)+i;
        D(piu)=1;
        DES=[DES D];
        end
    case {'bal'}
        %Base line vs others
        DES=-ones(1,C-1);
        DES=[DES;eye(C-1)];
    case {'trd'}
        for i=C-1:-1:1
            DES(1:i,i)=-1./length(1:i);
            DES(i+1:C,i)=1./length(i+1:C);
        end 
    otherwise
            disp('Unknown TYPE method. Resulting DES matrix is empty.')
            DES=[];
end
        