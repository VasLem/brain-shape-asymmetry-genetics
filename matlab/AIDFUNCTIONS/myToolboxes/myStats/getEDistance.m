function D = getEDistance(A,B)
         D = sqrt(sum((A-B).^2,2));
end