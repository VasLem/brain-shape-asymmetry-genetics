function [ind12,ind21,raw] = SNPPoslookup(POS1,POS2)   
   [~,ind12] = ismember(POS1(:),POS2(:));
   raw = ind12;
   ind21 = find(ind12);
   ind12 = ind12(ind21);
end