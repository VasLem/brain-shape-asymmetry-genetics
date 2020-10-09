function out = crop(obj,varargin)
%reduces 'mesh' so that it only contains the vertices indicated in 'indices'
%	: indices : vector containing the indices of the vertices you want to
%	keep or delete dependent on the action in varargin
    if nargout == 1
           obj = clone(obj);
           out = obj;
    end
    Input = find(strcmp(varargin, 'Action'));
    if isempty(Input)
        action = 'crop';
    else
        action = varargin{Input+1};
    end
    [Vindex,Findex] = getVindexFindex(obj,varargin{:});
    switch action
        case 'crop'
            %if length(Vindex)==size(obj.Vertices,2), return; end
        case 'delete'
            if isempty(Vindex), return; end
            fullindex = (1:1:size(obj.Vertices,2));
            Vindex = setdiff(fullindex,Vindex);
            Findex = Vindex2Findex(obj,Vindex);
    end
    
    %throw away all quads containing an index not in the indices we want to keep
    [tmp,LOC] = ismember(obj.Faces,Vindex);
    
    obj.Faces = LOC(:,Findex);
    obj.Vertices = obj.Vertices(:,Vindex);
    updateFields(obj,'Reduce',Vindex);
end
    

    
