function updateModel(obj)
            if isempty(obj.Model)
                ActivateAll(obj,'off');
                return;
            end
            % Connecting the faces
            obj.BaseCoeff = zeros(1,obj.Model.nrEV);
            obj.CurrentFace = clone(obj.Model.Average);
            obj.CurrentCoeff = zeros(1,obj.Model.nrEV);
            obj.CurrentFace.Axes = obj.Viewer.RenderAxes;
            obj.CurrentFace.Selected = true;
            obj.CurrentFace.Visible = true;
            % taking care of PC panel
            ActivatePCPanel(obj,'on');
            ActivateCARPanel(obj,'on');
            ActivateButton(obj,'on');
            for i=1:1:obj.Model.nrEV
                list{i} = num2str(i); %#ok<AGROW>
            end
            set(obj.Handles.PC.Popupmenu,'string',list);
            updatePCPanel(obj,1);
            updateCARPanel(obj);
            if ~strcmp(obj.Model.Type,'propertyPCA'), return; end
            ActivatePROPPanel(obj,'on');
            set(obj.Handles.PROP.Popupmenu,'String',obj.Model.PropNames);
            set(obj.Handles.PROP.Listbox,'String',obj.Model.PropNames);
            set(obj.Handles.PROP.Listbox,'max',length(obj.Model.PropNames));
            set(obj.Handles.PROP.Listbox,'min',0);
            set(obj.Handles.PROP.Listbox,'Value',[]);
            updatePROPPanel(obj,1);
end