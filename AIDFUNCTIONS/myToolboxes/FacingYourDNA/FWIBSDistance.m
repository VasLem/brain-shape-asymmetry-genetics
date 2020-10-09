function out = FWIBSDistance(v1,v2,fA,fa)
         ind = find(~isnan(v1));
         v1 = v1(ind);v2 = v2(:,ind);fA = fA(ind);fa = fa(ind);
         
         NfA = repmat(1./fA,s,1);
         NfA(isnan(v2)) = 0;
         Nfa = repmat(1./fa,s,1);
         Nfa(isnan(v2)) = 0;    
         N = sum(NfA+Nfa,2);
         
         
         
         N = 2*ones(size(v2));
         N(isnan(v2)) = 0;
         N = sum(N,2);
         
         A1 = v1;A1(v1==-1) = 2;A1(v1==0) = 1;A1(v1==1) = 0;%A1(isnan(v1)) = 0;
         A2 = v2;A2(v2==-1) = 2;A2(v2==0) = 1;A2(v2==1) = 0;%A2(isnan(v2)) = 0;
         
         IBSA = zeros(size(A2));
         A12 = A1==2;
         IBSA(:,A12) = (A2(:,A12)-2)+2;
         A11 = A1==1;
         tmp = (A2(:,A11)-1)+1;
         tmp(tmp==2) = 1;
         IBSA(:,A11)= tmp;
         
         A1 = v1;A1(v1==-1) = 0;A1(v1==0) = 1;A1(v1==1) = 2;%A1(isnan(v1)) = 0;
         A2 = v2;A2(v2==-1) = 0;A2(v2==0) = 1;A2(v2==1) = 2;%A2(isnan(v2)) = 0;
         
         IBSa = zeros(size(A2));
         A12 = A1==2;
         IBSa(:,A12) = (A2(:,A12)-2)+2;
         A11 = A1==1;
         tmp = (A2(:,A11)-1)+1;
         tmp(tmp==2) = 1;
         IBSa(:,A11)= tmp;
         
         out = 1-nansum((IBSA./repmat(fA,s,1))+(IBSa./repmat(fa,s,1)),2)./N;      
         
                  IBSA./repmat(fA,s,1)
         
         
end