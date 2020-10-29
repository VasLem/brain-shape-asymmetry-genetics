function [N,D] = KNN(obj,obj2,k)
% Finds the k nearest neighbors of obj2 to obj, uses kde toolbox
    switch class(obj2)
        case {'meshObj' 'LMObj'}
            pl = obj2.Vertices;
        otherwise
            pl = obj2;
    end
    [N,D] = knn(kde(obj.Vertices,5),pl,k);% kde toolbox routine
    if k>1
       D = zeros(k,length(pl));
       for i=1:1:k
           D(i,:) = sqrt(sum((obj.Vertices(:,N(i,:))-pl).^2));
       end
    end

end