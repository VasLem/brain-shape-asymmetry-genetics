function updateToolbar(obj,modetype,oldmode)
    switch modetype
        case 'Mode'
            switch oldmode
                case 'none'
                    return;
                case 'light'
                    if ~strcmp(obj.Mode,'light'), set(obj.Toolbar.light_mode,'State','off'); end
                case 'camera'
                    if ~strcmp(obj.Mode,'camera'), set(obj.Toolbar.cam_mode,'State','off'); end
                otherwise
                    return;
            end
        case 'SelectionMode'
            switch oldmode
                case 'none'
                    return;
                case 'landmark'
                    if ~strcmp(obj.SelectionMode,'landmark'), set(obj.Toolbar.landmark_mode,'State','off'); end
                case 'area'
                    if ~strcmp(obj.SelectionMode,'area'), set(obj.Toolbar.area_mode,'State','off'); end
                case 'fill'
                    if ~strcmp(obj.SelectionMode,'fill'), set(obj.Toolbar.fill_mode,'State','off'); end
                case 'brush'
                    if ~strcmp(obj.SelectionMode,'brush'), set(obj.Toolbar.brush_mode,'State','off'); end
                case 'full'
                    if ~strcmp(obj.SelectionMode,'full'), set(obj.Toolbar.full_mode,'State','off'); end
                otherwise
                    return;
            end
        otherwise
            return;
    end
end