function val = getLabelResultsIntoVoxelValuesv2(results,LABELS)
   val = nan*zeros(size(LABELS),'single')';
   nLC = length(results);
   HI = HierarchicalInterface;
   HI.nL = size(LABELS,2);
   for i=1:1:nLC
       [L,C] = Ind2LC(HI,i);
       lab = LABELS(:,L)';
       val(L,lab==C) = results(i);
   end
end