function updateVoxelImage(obj,action,key,varargin)
    if ~isempty(find(strcmp(varargin,'Busy'),1)),obj.Status = 'Busy';drawnow;end
    children = obj.VoxelImageChildren;
    for i=1:1:length(children)
        switch action
            case 'Update SliceMode'
                mode = children{i}.SliceMode;
                switch key
                    case 'x'
                        if ~isempty(strfind(mode,'x'))
                           mode = setdiff(mode,'x');
                        else
                           mode = union(mode,'x');
                        end
                    case 'y'
                        if ~isempty(strfind(mode,'y'))
                           mode = setdiff(mode,'y');
                        else
                           mode = union(mode,'y');
                        end
                    case 'z'
                        if ~isempty(strfind(mode,'z'))
                           mode = setdiff(mode,'z');
                        else
                           mode = union(mode,'z');
                        end
                    otherwise
                        continue;
                end
                children{i}.SliceMode = mode;
            case 'Update SliceNr'
                switch key
                    case 'uparrow'
                        val = +1;
                    case 'downarrow'
                        val = -1;
                    otherwise
                        continue;
                end
                modifier = varargin{1};
                if isempty(modifier)
                   children{i}.ZSliceNr = children{i}.ZSliceNr +val;
                elseif strcmp(modifier{1},'alt')
                   children{i}.XSliceNr = children{i}.XSliceNr +val; 
                elseif strcmp(modifier{1},'shift')
                   children{i}.YSliceNr = children{i}.YSliceNr +val; 
                else
                    continue;
                end                
            otherwise
        end
    end
    obj.Status = 'Ready';
    if isempty(obj.Parent), return; end
    handles = guidata(obj.Parent);
    handles.viewerChange(handles,'VoxelImage Updated');
end