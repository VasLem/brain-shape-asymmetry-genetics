function importStruct(obj,filename,varargin)
    % changing directory
         Input = find(strcmp(varargin,'Path'));
         if ~isempty(Input)
            cd(varargin{Input+1});
         end
     % Getting type
         Input = find(strcmp(varargin,'Type'));
         if isempty(Input)
            type = 'Mat Struct';
         else
            type = varargin{Input+1};
         end
     % loading    
       in = load(filename);
       loaded = fields(in);
       %if strcmp(class(in.(loaded{1})),'meshObj');
       struct = in.(loaded{1});
       struc2obj(obj,struct);
%        copy(in.(loaded{1}),obj);
%        delete(in.(loaded{1}));
end