function varargout = viewer(obj,varargin)

    Input = find(strcmp(varargin, 'Viewer'));
    if ~isempty(Input)
        if strcmp(class(varargin{Input+1}),'viewer3DObj'), v = varargin{Input+1};end
    else
        Input = find(strcmp(varargin, 'CalledProcess'));
        if ~isempty(Input)
           v = viewer3DObj('Tag','3D Viewer','CalledProcess',varargin{Input+1});
        else
           v = viewer3DObj('Tag','3D Viewer');
        end
    end
    obj.Axes = v.RenderAxes;
    obj.Visible = true;
    obj.Selected = true;
    if nargout == 1
        varargout{1} = v;
    end
end