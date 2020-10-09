function CM = combineFMatches(M1,M2)
         CM = zeros(size(M1,1),size(M2,2),size(M1,3)+size(M2,3));
         CM(:,:,1:size(M1,3)) = M1;
         CM(:,:,size(M1,3)+1:end) = M2;
end