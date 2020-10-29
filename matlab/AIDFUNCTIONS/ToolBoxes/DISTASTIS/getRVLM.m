function RV = getRVLM(RF)
         %[nLM,nV,nO] = size(RF);
         
     nfaces = size(RF,1);
     nV = size(RF,2);
     njuges = size(RF,3);
     %S = RF;
     ncells=nV*njuges;
     Grosse_mat=zeros(ncells,nfaces);
     %f = waitbar(0,'RV matrix');
     for k=1:nfaces;
         %k=1;
         la_c = squeeze(RF(k,:,:));
         Grosse_mat(:,k)=la_c(:);
         %waitbar(k/njuges,f);
     end
     %close(f);
     num_RV=Grosse_mat'*Grosse_mat;
     la_diag=diag(num_RV);
     nd=length(la_diag);
     den_RV=(repmat(la_diag',nd,1).*repmat(la_diag,1,nd)).^(1/2);
     RV=num_RV./den_RV; 
end