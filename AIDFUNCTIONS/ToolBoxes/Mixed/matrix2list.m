function list = matrix2list(matrix,nr)
         list = reshape(matrix,nr,size(matrix,2)*size(matrix,3));
end