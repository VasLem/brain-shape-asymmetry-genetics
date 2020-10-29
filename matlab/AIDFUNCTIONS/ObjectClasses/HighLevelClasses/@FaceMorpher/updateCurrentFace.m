function updateCurrentFace(obj,what)
         switch what
             case 'PC'
                 pc = get(obj.Handles.PC.Popupmenu,'Value');
                 value = get(obj.Handles.PC.Slider,'Value');
                 obj.CurrentCoeff(pc) = value*obj.Model.EigStd(pc);
                 updatingFace(obj);
                 updatePROPPanel(obj);
                 updateCARPanel(obj);
             case 'PROP'
                 prop = get(obj.Handles.PROP.Popupmenu,'Value');
                 value = get(obj.Handles.PROP.Slider,'Value');
                 indprop = get(obj.Handles.PROP.Listbox,'Value');
                 indprop = setdiff(indprop,prop);
                 path = getPropertyPath(obj.Model,prop,indprop);
                 obj.CurrentCoeff = setPropertyValue(obj.Model,obj.CurrentCoeff,prop,value,path);
                 updatingFace(obj);
                 updatePCPanel(obj);
                 updateCARPanel(obj);
             case 'CAR'
                 value = get(obj.Handles.CAR.Slider,'Value');
                 if get(obj.Handles.CAR.Popupmenu,'Value')==2
                    value = -1*value;
                 end
                 obj.CurrentCoeff = value*obj.CarCoeff;
                 updatingFace(obj);
                 updatePCPanel(obj);
                 updatePROPPanel(obj);
             case 'RESET'
                 obj.CurrentCoeff = obj.BaseCoeff;
                 updatingFace(obj);
                 updatePCPanel(obj);
                 updatePROPPanel(obj);
                 updateCARPanel(obj);
             otherwise 
                 return;
         end
end

function updatingFace(obj)
         scan = getScan(obj.Model,obj.CurrentCoeff);
         obj.CurrentFace.Vertices = scan.Vertices;
         
         if isempty(scan.TextureMap)
             obj.CurrentFace.TextureColor = scan.TextureColor;
         else
             obj.CurrentFace.TextureMap = clone(scan.TextureMap);
         end
         obj.CurrentFace.UserData = scan.UserData;
         basescan = getScan(obj.Model,obj.BaseCoeff);
%          if ~isempty(scan.Value)
%             obj.CurrentFace.Value = scan.Value;
%          else
%            obj.CurrentFace.Value = vDistances(obj.CurrentFace,basescan);
%          end
         delete(scan);
         delete(basescan);
end