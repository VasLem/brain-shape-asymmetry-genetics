function out = getASDistance(v1,v2)
%          v1IBA = v1;v1IBA(v1==-1) = 2;v1IBA(v1==0) = 1;v1IBA(v1==1) = 0;
%          v1IBa = v1;v1IBa(v1==-1) = 0;v1IBa(v1==0) = 1;v1IBa(v1==1) = 2;
%          
%          v2IBA = v2;v2IBA(v2==-1) = 2;v2IBA(v2==0) = 1;v2IBA(v2==1) = 0;
%          v2IBa = v2;v2IBa(v2==-1) = 0;v2IBa(v2==0) = 1;v2IBa(v2==1) = 2;
%          
%          IBSA = (v1IBA+v1IBA2);
         out = nanmean(abs(repmat(v1,size(v2,1),1)-v2),2);
end