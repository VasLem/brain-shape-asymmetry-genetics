function [RV] = getRvCoeff(S)
     nfaces = size(S,1);
     njuges = size(S,3);
     ncells=nfaces*nfaces;
     Grosse_mat=zeros(ncells,njuges);
     f = waitbar(0,'RV matrix');
     for k=1:njuges;
         la_c = S(:,:,k);
         Grosse_mat(:,k)=la_c(:);
         waitbar(k/njuges,f);
     end
     close(f);
     num_RV=Grosse_mat'*Grosse_mat;
     la_diag=diag(num_RV);
     nd=length(la_diag);
     den_RV=(repmat(la_diag',nd,1).*repmat(la_diag,1,nd)).^(1/2);
     RV=num_RV./den_RV; 
end