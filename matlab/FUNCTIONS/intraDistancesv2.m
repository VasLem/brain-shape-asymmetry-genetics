        function out = intraDistancesv2(obj,varargin)
            [Vindex,~] = getVindexFindex(obj,varargin{:});
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
            out = nan*zeros(1,obj.nVertices);
            To_Do = 1:obj.nVertices;
            To_Do = setdiff(To_Do,Vindex);
            current = Vindex;
            out(current) = 0;  
            Neighbors = 1;
            counter = 0; 
            while ~isempty(Neighbors)&&(counter<obj.nVertices)
                 counter = counter + 1;
                 [current_index,Neighbors,dist] = find(A(current,:));
                 dist = dist(:)' + out(current(current_index));
                 [~,I,~] = unique(Neighbors,'first');
                 dist_first = dist(I);
                 [Neighbors,I,~] = unique(Neighbors,'last');
                 dist_last = dist(I);
                 dist = min([dist_first;dist_last]);
                 [Neighbors,~,J] = intersect(To_Do,Neighbors);
                 dist = dist(J);
                 out(Neighbors) = dist;
                 if isempty(find(out(Neighbors)<=max_dist,1)), break; end
                 current = Neighbors;
                 To_Do = setdiff(To_Do,Neighbors);
            end
            out(find(out>max_dist,1)) = nan;
        end