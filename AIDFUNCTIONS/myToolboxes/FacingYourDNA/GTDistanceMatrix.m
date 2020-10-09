function D = GTDistanceMatrix(v,fAA,faA,faa)
         s = size(v,1);
         D = nan*zeros(s,s);
         for i=1:1:s
             D(i,:) = GTDistance(v(i,:),v,fAA,faA,faa)';
         end
end