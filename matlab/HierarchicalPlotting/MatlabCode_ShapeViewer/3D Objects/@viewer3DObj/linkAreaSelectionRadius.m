function linkAreaSelectionRadius(obj)
    % case 1: does not exist so create one
    if isempty(obj.AreaSelectionRadius)||isstruct(obj.AreaSelectionRadius)
       obj.AreaSelectionRadius = areaObj('Visible',false,'Parent',obj.CurrentMesh);
       obj.AreaSelectionRadius.SingleColor = [0 0.7 0];
       obj.AreaSelectionRadius.ViewMode = 'Wireframe';       
%        areaBorder(obj.AreaSelectionRadius);
%        obj.AreaSelectionRadius.Border.SingleColor = [0 1 0];
       return;
    end
    % case 2: existed before, but was deleted
    if ~validChild(obj,obj.AreaSelectionRadius)
       obj.AreaSelectionRadius = areaObj('Visible',false,'Parent',obj.CurrentMesh);
       obj.AreaSelectionRadius.SingleColor = [0 0.7 0];
       obj.AreaSelectionRadius.ViewMode = 'Wireframe';
%        areaBorder(obj.AreaSelectionRadius);
%        obj.AreaSelectionRadius.Border.SingleColor = [0 1 0];
       return;
    end
    % case 3: does exist but is not linked to any parent
    if isempty(obj.AreaSelectionRadius.Parent)
       obj.AreaSelectionRadius.Parent = obj.CurrentMesh;
       return;
    end
    % case 4: does exist but is not properly linked to the current mesh
    if ~(obj.AreaSelectionRadius.Parent==obj.CurrentMesh)
       delete(obj.AreaSelectionRadius);
       obj.AreaSelectionRadius = areaObj('Visible',false,'Parent',obj.CurrentMesh);
       obj.AreaSelectionRadius.SingleColor = [0 0.7 0];
       obj.AreaSelectionRadius.ViewMode = 'Wireframe';
%        areaBorder(obj.AreaSelectionRadius);
%        obj.AreaSelectionRadius.Border.SingleColor = [0 1 0];
       return;
    end
end