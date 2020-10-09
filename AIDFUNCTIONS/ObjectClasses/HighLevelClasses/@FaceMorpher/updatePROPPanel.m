function updatePROPPanel(obj,prop)
      try
         if nargin < 2, 
            prop = get(obj.Handles.PROP.Popupmenu,'Value');
         else
            set(obj.Handles.PROP.Popupmenu,'Value',prop);
         end
         value = obj.CurrentFace.UserData.Prop(prop);
         set(obj.Handles.PROP.Edit,'String',num2str(value));
         propvalues = repmat(obj.Model.AvgProp(prop),1,obj.Model.n) + obj.Model.EigVec(end-obj.Model.nrP+prop,:)*obj.Model.Tcoeff';
         spread = std(propvalues);
         set(obj.Handles.PROP.Slider,'Min',min(propvalues)-obj.PropScale*spread,'Max',max(propvalues)+obj.PropScale*spread);       
         set(obj.Handles.PROP.Slider,'Value',value);
      catch
      end
end