function varargout = figure(obj,varargin)
    if isempty(varargin)
       f = figure;axis equal;
    elseif ~strcmp(get(varargin{1},'Type'),'figure')
       f = figure;axis equal; 
    else
       f = varargin{1};
       figure(f);axis equal;
    end
    obj.Axes = gca;
    obj.Visible = true;
    if nargout == 1
        varargout{1} = f;
    end
end