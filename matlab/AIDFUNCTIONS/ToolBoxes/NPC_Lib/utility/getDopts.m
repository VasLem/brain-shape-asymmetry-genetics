function [D, options]=getDopts(D,dim,options)
% return the matrix data D and labels. D in iNPut char, cell or struct 

if nargin<3, 
	options=1; 
else
    if isnumeric(options), temp=options; clear options; options.OUT=temp; 
    elseif  isstruct(options) if not(isfield(options,'OUT')) options.OUT=1; end 
	end
end

if isnumeric(options), temp=options; clear options; options.OUT=temp; clear temp; end

if ischar(D) 
    [no D]=v({D});
    options.labels.dims{dim}=[D.name];
    options.labels.dimslabel{dim}='';
    [D]=[D.vals];
elseif iscell(D)
    [no D]=v(D);
    options.labels.dims{dim}=[D.name];
    options.labels.dimslabel{dim}='';
    [D]=[D.vals];
elseif isstruct(D)
    options.labels.dims{dim}=[D.name];
    options.labels.dimslabel{dim}='';
    [D]=[D.vals];
end

if isvector(D)
    D=D(:);    
end
