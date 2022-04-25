function exportStruct(obj,filename,varargin) %#ok<INUSL>
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Opening file and write Heading
         struct = obj2struc(obj);
         save(filename,'struct');
end