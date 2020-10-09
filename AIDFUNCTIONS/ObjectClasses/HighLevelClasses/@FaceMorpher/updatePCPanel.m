function updatePCPanel(obj,pc)
         if nargin < 2, 
            pc = get(obj.Handles.PC.Popupmenu,'Value');
         else
            set(obj.Handles.PC.Popupmenu,'Value',pc);
         end
         value = obj.CurrentCoeff(pc)/obj.Model.EigStd(pc);
         set(obj.Handles.PC.Edit,'String',num2str(value));
         set(obj.Handles.PC.Slider,'Value',value);
end