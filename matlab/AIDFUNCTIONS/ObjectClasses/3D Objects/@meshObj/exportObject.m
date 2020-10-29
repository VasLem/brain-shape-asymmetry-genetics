function exportObject(obj,filename,varargin) %#ok<INUSL>
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Opening file and write Heading
         save(filename,'obj');
end