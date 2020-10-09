function out = checkVerticesColorUpdate(obj,action)
out = false;
switch action
    case 'update texture'
        if strcmp(obj.ColorMode,'Texture')&& (size(obj.TextureColor,2)==size(obj.Vertices,2)), setPatch(obj,'VerticesColor'); out = true; end
    case 'update index'
        if strcmp(obj.ColorMode,'Indexed')&& (length(obj.IndexedColor)==size(obj.Vertices,2)), setPatch(obj,'VerticesColor'); out = true;end
    case 'update vertices'
        switch obj.ColorMode
               case 'Single'
                    setPatch(obj,'VerticesColor'); out = true;          
               case 'Indexed'
                    if (size(obj.Vertices,2)==length(obj.IndexedColor)), setPatch(obj,'VerticesColor'); out = true; end
               case 'Texture'
                    if (size(obj.Vertices,2)==size(obj.TextureColor,2)), setPatch(obj,'VerticesColor'); out = true; end
        end
    otherwise
end
end