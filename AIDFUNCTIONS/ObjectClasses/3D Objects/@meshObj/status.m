function status = status(obj)
            status{1} = obj.Selected;
            status{2} = obj.Visible;
            status{3} = obj.Tag;
            if validChild(obj,obj.PoseLM)
                status{4} = true;
            else
                status{4} = false;
            end
            status{5} = obj.Mapped;
end