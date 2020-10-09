function out = getMatFile(string)
         if nargin < 1,string = 'select file';end
         [filename, pathname] = uigetfile({'*.mat','Mat file'},string,'MultiSelect','on');
         if isequal([filename,pathname],[0,0]), out = []; return; end
         cd(pathname);
         if iscell(filename)
                out = cell(1,length(filename));
                for i=1:1:length(filename)
                    in = load(filename{i});
                    loaded = fields(in);
                    out{i} = in.(loaded{1});
                end
         else
             out = cell(1,1);
             in = load(filename);
             loaded = fields(in);
             out{1} = in.(loaded{1});
         end
end