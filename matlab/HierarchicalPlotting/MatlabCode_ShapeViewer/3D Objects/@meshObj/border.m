function varargout = border(obj,varargin)
% returns the border points of a full mesh or a selection of vertices and
% faces
try
       [Vindex,Findex] = getVindexFindex(obj,varargin{:});
       clear Findex;
       Area = areaObj('Visible',false,'Parent',obj,'VerticesIndex',Vindex);
       B = border(Area);
       if isempty(B), varargout{1} = []; return; end
       if nargout == 1
          varargout{1} = B;
       else
          obj.Border = B;
          obj.Border.VerticesIndex = Area.VerticesIndex(obj.Border.VerticesIndex);
          updateChildren(obj,'Selected Change');
          obj.Border.Axes = obj.Axes;
          obj.Border.Visible = obj.Visible;          
       end
       delete(Area);
catch
    disp('Failed to detect the border');
end
end
