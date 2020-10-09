function D = pairWiseGeodesicDistancesFromAdjacency(A)
    nrV = size(A,1);
    D = zeros(nrV,nrV);
    [path,ID] = setupParForProgress(nrV);
    parfor k=1:nrV
        D(k,:) = subfunction(A,k);
        parfor_progress;
    end
    closeParForProgress(path,ID);
end
function [distances] = subfunction(A,Vindex)
    nrV = size(A,1);
    max_dist = inf;
    %initialisatie
    distances = nan*zeros(1,nrV);
    To_Do = (1:nrV);
    To_Do = setdiff(To_Do,Vindex);
    current = Vindex;
    distances(current) = 0;  
    Neighbors = 1;
    counter = 0; 
    while ~isempty(Neighbors)&&(counter<nrV);
         counter = counter + 1;
         [current_index,Neighbors,dist] = find(A(current,:));
         dist = dist(:)' + distances(current(current_index));
         [~, I, ~] = unique(Neighbors,'first');
         dist_first = dist(I);
         [Neighbors, I, ~] = unique(Neighbors,'last');
         dist_last = dist(I);
         dist = min([dist_first;dist_last]);
         [Neighbors, ~, J] = intersect(To_Do,Neighbors);
         dist = dist(J);
         distances(Neighbors) = dist;
         if isempty(find(distances(Neighbors)<=max_dist,1)), break; end
         current = Neighbors;
         To_Do = setdiff(To_Do,Neighbors);
    end
    distances(find(distances>max_dist,1)) = nan;
end