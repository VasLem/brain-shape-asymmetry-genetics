function [distances] = intraDistances(obj,varargin)
    [Vindex,Findex] = getVindexFindex(obj,varargin{:});
    Input = find(strcmp(varargin, 'Type'));
    if ~isempty(Input)
       distance_type = varargin{Input+1};
    else
       distance_type = 'triangle';
    end
    Input = find(strcmp(varargin, 'Max'));
    if ~isempty(Input)
       max_dist = varargin{Input+1};
    else
       max_dist = inf;
    end

    switch distance_type
        case 'triangle'
            A = obj.Adjacency;
        case 'edge'
            A = vertexAdjacency(obj);
        otherwise
            return;
    end
    %initialisatie
        distances = nan*zeros(1,size(obj.Location,2));
        To_Do = (1:size(obj.Location,2));
        To_Do = setdiff(To_Do,Vindex);
        current = Vindex;
        distances(current) = 0;  
        Neighbors = 1;
        counter = 0; 
    while ~isempty(Neighbors)&&(counter<size(obj.Location,2));
         counter = counter + 1;
         [current_index,Neighbors,dist] = find(A(current,:));
         dist = dist(:)' + distances(current(current_index));
         [tmp, I, J] = unique(Neighbors,'first');
         dist_first = dist(I);
         [Neighbors, I, J] = unique(Neighbors,'last');
         dist_last = dist(I);
         dist = min([dist_first;dist_last]);
         [Neighbors, I, J] = intersect(To_Do,Neighbors);
         dist = dist(J);
         distances(Neighbors) = dist;
         if isempty(find(distances(Neighbors)<=max_dist,1)), break; end
         current = Neighbors;
         To_Do = setdiff(To_Do,Neighbors);
    end
    distances(find(distances>max_dist,1)) = nan;
end
