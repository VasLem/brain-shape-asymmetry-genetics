function vargout = getPatch(obj,varargin)
         if isempty(obj.ph) || ~ishandle(obj.ph), vargout{1} = []; return; end
         vargout = get(obj,ph,varargin{:});
end % get Patch Properties.