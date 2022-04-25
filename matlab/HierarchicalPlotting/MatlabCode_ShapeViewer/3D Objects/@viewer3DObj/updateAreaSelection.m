function updateAreaSelection(obj,action,varargin)
        switch action
            case 'Add Area'
                if isempty(varargin{1}),return; end
                obj.AreaSelection.VerticesIndex = union(obj.AreaSelection.VerticesIndex, varargin{1});
            case 'Remove Area'
                if isempty(varargin{1}),return; end
                obj.AreaSelection.VerticesIndex = setdiff(obj.AreaSelection.VerticesIndex, varargin{1});
            case 'Intersect Area'
                if isempty(varargin{1}),return; end
                obj.AreaSelection.VerticesIndex = intersect(obj.AreaSelection.VerticesIndex, varargin{1});
            case 'Set Area'
                if isempty(varargin{1}),return; end
                obj.AreaSelection.VerticesIndex = varargin{1};
            case 'Delete'
                if ~isempty(obj.CurrentMesh)
                    crop(obj.CurrentMesh,'VertexIndex',obj.AreaSelection.VerticesIndex,'Action','delete');
                end
                obj.AreaSelection.Visible = false;
                obj.AreaSelection.VerticesIndex = [];
            case 'Crop'
                if ~isempty(obj.CurrentMesh)
                   crop(obj.CurrentMesh,'VertexIndex',obj.AreaSelection.VerticesIndex,'Action','crop');
                end
                obj.AreaSelection.Visible = false;
                obj.AreaSelection.VerticesIndex = [];
            case 'Clear'
                obj.AreaSelection.Visible = false;
                obj.AreaSelection.VerticesIndex = [];
            case 'Invert'
                if ~isempty(obj.CurrentMesh)
                   fullindex = (1:1:size(obj.CurrentMesh.Location,2));
                   obj.AreaSelection.VerticesIndex = setdiff(fullindex,obj.AreaSelection.VerticesIndex);
                end
            case 'Smooth Surface'
                if isempty(obj.AreaSelection.VerticesIndex), return; end
                obj.Status = 'Busy';drawnow;
                smoothSurface(obj.AreaSelection,obj.SmoothRuns,obj.SmoothMode);
                setPatch(obj.AreaSelection,'VerticesColor');
                obj.Status = 'Ready';
            case 'Smooth Color'
                if isempty(obj.AreaSelection.VerticesIndex), return; end
                obj.Status = 'Busy';drawnow;
                smoothColor(obj.AreaSelection,obj.SmoothRuns,obj.SmoothMode);
                setPatch(obj.AreaSelection,'VerticesColor');
                obj.Status = 'Ready';
            case 'Subdivide Triangles'
                if isempty(obj.AreaSelection.VerticesIndex), return; end
                subdivideTriangles(obj.AreaSelection,obj.SubdivideMode,obj.SubdivideVal);
            case 'Reduce Triangles'
                if isempty(obj.AreaSelection.VerticesIndex), return; end
                reduceTriangles(obj.AreaSelection,(100-obj.ReducePercentage)/100);
                obj.AreaSelection
            otherwise
        end
        if ~isempty(obj.AreaSelection)
             if length(obj.AreaSelection.VerticesIndex) > 2
                obj.AreaSelection.Visible = true;
                render(obj.AreaSelection);
             else
                obj.AreaSelection.Visible = false;
             end
        end
end