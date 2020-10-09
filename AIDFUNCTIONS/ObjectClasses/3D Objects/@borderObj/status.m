function status = status(obj)
            status{1} = obj.Selected;
            status{2} = obj.Visible;
            if ~isempty(obj.Parent)
                status{3} = obj.Parent.Tag;
            else
                status{3} = '';
            end
            status{4} = obj.Tag;
end