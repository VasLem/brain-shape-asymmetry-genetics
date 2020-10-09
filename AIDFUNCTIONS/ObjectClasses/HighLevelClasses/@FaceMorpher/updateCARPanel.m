function updateCARPanel(obj)
         obj.CarCoeff = obj.CurrentCoeff;
         switch get(obj.Handles.CAR.Popupmenu,'Value')
             case 1
                 val = 1;
             case 2
                 val = 0;
         end                
         set(obj.Handles.CAR.Edit,'String',num2str(val));
         set(obj.Handles.CAR.Slider,'Value',val);
end