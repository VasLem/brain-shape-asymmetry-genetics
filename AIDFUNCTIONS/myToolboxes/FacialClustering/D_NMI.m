function NMI = D_NMI(clustering1, clustering2)
% This function returns the normalised mutual information between two clusterings
 
    nClusters1 = max(clustering1);
    nClusters2 = max(clustering2);
    
    % create the contigency table
    CTable = zeros(nClusters1, nClusters2);
    for i = 1:nClusters1,
        ind_cl1 = double(clustering1==i);
        for j = 1:nClusters2,
            ind_cl2 = double(clustering2==j);
            CTable(i,j) = ind_cl1' * ind_cl2;
        end
    end
    n1 = sum(CTable,2);
    n2 = sum(CTable,1);
    n = sum(n1);   % or sum(n2);
    
    % NMI
    teller = 0;
    noemer1 = 0;
    noemer2 = 0;
    
    for i = 1:nClusters1,
        for j = 1:nClusters2,
            if CTable(i,j) ~= 0,
                teller = teller + (CTable(i,j) * log((n * CTable(i,j))/(n1(i)*n2(j))));
            end
            if i == 1,
                noemer2 = noemer2 + (n2(j) * log(n2(j)/n));
            end
        end
        noemer1 = noemer1 + (n1(i) * log(n1(i)/n));
    end 
    
    NMI = teller/sqrt(noemer1 * noemer2);
    
end

