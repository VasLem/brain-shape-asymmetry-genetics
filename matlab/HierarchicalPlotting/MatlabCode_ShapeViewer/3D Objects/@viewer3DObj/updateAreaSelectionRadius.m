function updateAreaSelectionRadius(obj, action, varargin)
        switch action
            case 'Add Area'
                if isempty(varargin{1}), updateAreaSelectionRadius(obj,'Clear'), return; end
                obj.AreaSelectionRadius.VerticesIndex = varargin{1};
            case 'Clear'
                obj.AreaSelectionRadius.Visible = false;
                obj.AreaSelectionRadius.VerticesIndex = [];
            otherwise
        end
        if ~isempty(obj.AreaSelectionRadius)&&isobject(obj.AreaSelectionRadius)
             if size(obj.AreaSelectionRadius.VerticesIndex,2) > 2
                obj.AreaSelectionRadius.Visible = true;
                render(obj.AreaSelectionRadius);
             end
        end
end