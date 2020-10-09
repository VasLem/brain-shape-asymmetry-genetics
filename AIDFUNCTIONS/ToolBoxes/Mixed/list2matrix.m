function matrix = list2matrix(list,nr,nrP)
         matrix = reshape(list,nr,size(list,2)/nrP,nrP);
end