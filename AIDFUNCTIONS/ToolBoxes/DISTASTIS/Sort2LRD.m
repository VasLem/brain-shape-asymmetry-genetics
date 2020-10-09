function [L,R,D] = Sort2LRD(sort)
% convert sorting results into Indicator matrices, Co-Occurence matrices an
% d distance matrices.
         [nobjects,njuges]=size(sort);
         L = cell(1,njuges); % Indicator matrix
         D = zeros(nobjects,nobjects,njuges);% Distance matrix
         R = D;% co-occurence matrix
         for k=1:njuges;
             col= sort(:,k);% getting the relevant sorting results for current judge
             ncol=max(col);% determine the amount of groups 
             Lwork=zeros(nobjects,ncol);% initialize
             work=zeros(nobjects,nobjects);
             for j=1:nobjects;
                 Lwork(j,col(j))=1;
             end
             L{k} = Lwork;
             for j1=1:nobjects-1;
                 for j2=j1+1:nobjects;
                     work(j1,j2)=(col(j1)~=col(j2));
                 end
             end
             work=work+work';
             D(:,:,k) = work;
             R(:,:,k) = 1-work;
         end
end