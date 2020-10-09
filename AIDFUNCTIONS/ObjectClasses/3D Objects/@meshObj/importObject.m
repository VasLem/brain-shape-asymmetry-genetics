function importObject(obj,filename,path)
       if nargin<3, path = pwd; end
       if ~strcmp(path(end),'/'), path = [path '/'];end
       filename = [path filename];   
       in = load(filename);
       loaded = fields(in);
       copy(in.(loaded{1}),obj);
       delete(in.(loaded{1}));
end