function varargout = border(obj)
% TO DO: first check on floating points
        vertex_index = (1:1:size(obj.Location,2));
        %triangle_index = unique(obj.Faces(:));
        %index = setdiff(vertex_index,triangle_index);
        N = hist(obj.Tri(:),vertex_index); % nr triangles per vertex;
        E = full(sum(obj.Adjacency)); % nr of adjecant elements per vertex
        %index = find(N==0);
        index = find(N<E);
        if isempty(index), varargout{1} = []; return; end
        Border = borderObj('Parent',obj,'VerticesIndex',index,'Visible',obj.Visible);
        if nargout == 1
           varargout{1} = Border;
        else
            obj.Border = Border;
            %if isempty(obj.ph)||~ishandle(obj.ph)||~obj.UpdatePatch, return; end
            %createPatch(obj.Border);
        end
end
