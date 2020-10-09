function updateChildren(obj,action,varargin)

switch action
    case 'Axe Change'
        if obj.LinkChildren
                if ~isempty(obj.PoseLM), obj.PoseLM.Axes = obj.Axes; end
                if validChild(obj,obj.Border), obj.Border.Axes = obj.Axes; end
                if ~isempty(obj.CustomLM)
                    for l=1:1:length(obj.CustomLM)
                        if validChild(obj,obj.CustomLM{l})
                           obj.CustomLM{l}.Axes = obj.Axes;
                        end
                    end
                end
                if ~isempty(obj.Area)
                    for a=1:1:length(obj.Area)
                        if validChild(obj,obj.Area{a})
                           obj.Area{a}.Axes = obj.Axes;
                        end
                    end
                end
        end
    case 'Visible Change'
        if obj.LinkChildren
                if validChild(obj,obj.PoseLM), obj.PoseLM.Visible = obj.Visible; end
                if validChild(obj,obj.Border), obj.Border.Visible = obj.Visible; end
                if ~isempty(obj.CustomLM)
                    for l=1:1:length(obj.CustomLM)
                        if validChild(obj,obj.CustomLM{l})
                           obj.CustomLM{l}.Visible = obj.Visible;
                        end
                    end
                end
                if ~isempty(obj.Area)
                    for a=1:1:length(obj.Area)
                        if validChild(obj,obj.Area{a})
                           obj.Area{a}.Visible = obj.Visible;
                        end
                    end
                end
        end
    case 'Selected Change'
        if ~validChild(obj,obj.Border), return; end
        switch obj.Selected
        case true
           obj.Border.SingleColor = [0.5 0 0];
        case false
           obj.Border.SingleColor = [0.5 0.5 0.5];
        end
    case 'Patch Created'
        %try
            if validChild(obj,obj.PoseLM), refreshPatch(obj.PoseLM); end
            if ~isempty(obj.Border), refreshPatch(obj.Border); end
            if ~(isempty(obj.CustomLM))
                for i=1:1:length(obj.CustomLM)
                    if validChild(obj,obj.CustomLM{i}) , refreshPatch(obj.CustomLM{i}); end
                end
            end
            if ~(isempty(obj.Area))
                for i=1:1:length(obj.Area)
                    if validChild(obj,obj.Area{i}), refreshPatch(obj.Area); end
                end
            end
        %catch
        %end
    case 'Delete'
        if ~isempty(obj.PoseLM), delete(obj.PoseLM); end
        if ~isempty(obj.TextureMap), delete(obj.TextureMap); end
        if ~isempty(obj.Border), delete(obj.Border); end
        if ~(isempty(obj.CustomLM))
            for i=1:1:length(obj.CustomLM)
                if validChild(obj,obj.CustomLM{i}) , delete(obj.CustomLM{i}); end
            end
        end
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i}), delete(obj.Area); end
            end
        end
        if ~isempty(obj.RBF), delete(obj.RBF); end
        if ~isempty(obj.KDE), delete(obj.KDE); end
    case {'Crop' 'Reduce'}
        if ~isempty(obj.Border)
            delete(obj.Border);
            try 
                border(obj);
            catch
            end
        end
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i})
                   obj.Area.VerticesIndex = intersect(obj.Area.VerticesIndex,varargin{1}); 
                end
            end
        end
    case 'Smooth'
        if validChild(obj,obj.Border), setPatch(obj.Border,'VerticesColor'); end
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i})
                   setPatch(obj.Area{i},'VerticesColor');
                end
            end
        end
    case 'Subdivide Triangles'
        if validChild(obj,obj.Border), delete(obj.Border); border(obj); end
        % Area children will be lost!!!
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i})
                   delete(obj.Area{i});
                end
            end
        end
    case 'Append'
        if validChild(obj,obj.Border), delete(obj.Border); border(obj); end
        %if validChild(obj,obj.TextureMap), delete(obj.TextureMap); end
        % Area children will be lost!!!
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i})
                   delete(obj.Area{i});
                end
            end
        end
    case {'Transform Rigid' 'Transform Affine'}
        if validChild(obj,obj.PoseLM), transformRigid(obj.PoseLM,varargin{1}); end
        if validChild(obj,obj.Border), setPatch(obj.Border,'VerticesColor'); end
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i})
                   setPatch(obj.Area{i},'VerticesColor');
                end
            end
        end
    case {'Transform RBF'}
        if validChild(obj,obj.PoseLM), transformRBF(obj.PoseLM,varargin{1}); end
        if validChild(obj,obj.Border), setPatch(obj.Border,'VerticesColor'); end
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i})
                   setPatch(obj.Area{i},'VerticesColor');
                end
            end
        end
    case 'Transform' % object oriented transformation model
        if validChild(obj,obj.PoseLM), transform(varargin{1},obj.PoseLM); end
        if ~isempty(obj.Border), setPatch(obj.Border,'VerticesColor'); end
        if ~(isempty(obj.Area))
            for i=1:1:length(obj.Area)
                if validChild(obj,obj.Area{i})
                   setPatch(obj.Area{i},'VerticesColor');
                end
            end
        end
    otherwise
end
        



end