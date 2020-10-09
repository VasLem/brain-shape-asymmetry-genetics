function [ind12,ind21,multiple] = SNPAllPoslookup(POS1,POS2)  
   ind21 = [];
   ind12 = [];
   nPOS1 = length(POS1);
   multiple = zeros(1,nPOS1);
   for i=1:nPOS1
       %i=1;
       val = POS2-POS1(i);
       ind = find(val==0);
       if isempty(ind), continue;end
       N = length(ind);
       ind21 = [ind21 i*ones(1,N)];
       ind12 = [ind12 ind];
       if N<2, continue;end
       multiple(i) = true;
   end
end