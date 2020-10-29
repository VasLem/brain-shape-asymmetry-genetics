function setMatFile(string,var)
         if nargin < 1,string = 'select file';end
         [filename, pathname] = uiputfile({'*.mat','Mat file'},'Save As',string);
         if isequal([filename,pathname],[0,0]), return; end
         cd(pathname);
         save(filename,'var');
end