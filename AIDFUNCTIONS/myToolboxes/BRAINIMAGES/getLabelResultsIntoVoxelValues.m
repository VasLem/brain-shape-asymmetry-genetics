function val = getLabelResultsIntoVoxelValues(results,LABELS)
   val = nan*zeros(size(LABELS),'single')';
   nL = length(results);
   for l=1:1:nL
      lab = LABELS(:,l)';
      res = results{l};
      nC = length(res);
      for c=1:1:nC
         val(l,lab==c) = res(c);
      end
   end
end