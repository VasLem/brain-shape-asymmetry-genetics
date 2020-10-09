function removeActiveFilesAfterShutDown(varargin)
         Input = find(strcmpi(varargin, 'path'));
         if isempty(Input),mainpath = [pwd '/'];else mainpath = varargin{Input+1};end
         Input = find(strcmpi(varargin, 'id'));
         if isempty(Input),id = getComputerID;else id = varargin{Input+1};end
         allpaths = genpath(mainpath);
         ind = strfind(allpaths,':');
         ind = [1 ind];
         for i=2:1:length(ind)
            if i==2
               subpath = mainpath;
            else
               subpath = allpaths(ind(i-1)+1:ind(i)-1);
            end
            removeActiveFilesv2('path',subpath,'id',id);  
         end
end