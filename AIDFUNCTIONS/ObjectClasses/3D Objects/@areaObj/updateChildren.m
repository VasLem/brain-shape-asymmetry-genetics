function updateChildren(obj,action)

switch action
    case 'Axe Change'
        %if obj.LinkChildren
           if validChild(obj,obj.Border), obj.Border.Axes = obj.Axes; end
        %end
    case 'Visible Change'
        %if obj.LinkChildren
           if validChild(obj,obj.Border), obj.Border.Visible = obj.Visible; end
        %end
    case 'Index Change'
        if validChild(obj,obj.Border), delete(obj.Border); border(obj); end
    case 'Delete'
        if validChild(obj,obj.Border), delete(obj.Border); end
    otherwise
end
        



end