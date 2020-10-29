function setTitle(obj,String)
      if isempty(obj.Title)
         title(obj.RenderAxes,String,'Color',[1 1 1]);
      else
         set(obj.t,'String',String);
      end
end