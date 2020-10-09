function options=get_options(X,outputype,options)

%global NPCD__
if nargin==1, 
	outputype='NPC'; 
end
if nargin<3, 
	options=[]; 
end

if nargin<0
    X=1; 
elseif iscell(X)
    for i=1:length(X)
        dims(i)=max(1,length(X{i}));
    end
 
    options.labels=label_std(outputype,ones(dims));
    for i=1:length(X)
        options.labels.dims{i}=X{i};
    end
       X=ones(dims);
end



if nargin<3
    options.OUT=[1 0];
    options.Pobs=0;
    options.tail=0;
    options.Combdims=length(size(X));
else
    if isempty(options) options.OUT=1; end
    if isnumeric(options)
     temp=options; clear options; options.OUT=temp; clear temp;
     options.Pobs=0;
     options.tail=0;
     options.Combdims=length(size(X));

    else
        if isfield(options,'OUT')==0,     options.OUT=[1 0]; end
        if isfield(options,'Pobs')==0,    options.Pobs=0; end
        if isfield(options,'tail')==0,     options.tail=0; end
        if isfield(options,'Combdims')==0,     options.Combdims=length(size(X)); 
        %else options.Combdims=min(length(size(X)),options.Combdims);
        end
    end
end

if strmatch(outputype,'BF','exact')
        if isfield(options,'center')==0,     options.center='median'; end
end

if [strmatch(outputype,'ReM','exact')  strmatch(outputype,'MC','exact')]
    if isfield(options,'DES')==0,     options.DES=DES_ReM_std(size(X,1),'All'); 
    elseif ischar(options.DES),       options.DES=DES_ReM_std(size(X,1),options.DES); 
    elseif isempty(options.DES),      options.DES=DES_ReM_std(size(X,1),'All'); 
    end
    if isfield(options,'labels'), 
	templab=options.labels; 
	else
	templab=[]; 
	end  
        options.labels=label_std(outputype,X,options.DES);   
    if isfield(options,'labels'),
        if isfield(templab,'dims')
        for i=1:length(templab.dims)
            if not(isempty(templab.dims{i}))
                options.labels.dims{i}=templab.dims{i};
            end
        end
        end
        if isfield(templab,'dimslabel'),
        for i=1:length(templab.dimslabel)
            if not(isempty(templab.dimslabel{i}))
                options.labels.dimslabel{i}=templab.dimslabel{i};
            end
        end
        end
    end  
end


if isempty(options)
    options.labels=label_std(outputype,X);
else
    if isnumeric(options)
     options.labels=label_std(outputype,X);
    else
    if isfield(options,'labels'), 
	templab=options.labels; 
	else
	templab=[]; 
	end  
        options.labels=label_std(outputype,X,options);   
    if isfield(options,'labels'),
        if isfield(templab,'dims')
        for i=1:length(templab.dims)
            if not(isempty(templab.dims{i}))
                options.labels.dims{i}=templab.dims{i};
            end
        end
        end
        if isfield(templab,'dimslabel'),
        for i=1:length(templab.dimslabel)
            if not(isempty(templab.dimslabel{i}))
                options.labels.dimslabel{i}=templab.dimslabel{i};
            end
        end
        end
    end  

    end
end



if [ strmatch(outputype,'1s','exact') strmatch(outputype,'ReM','exact')]
    options.MU=0;
end

if strmatch(outputype,'StOrd','exact')
        if isfield(options,'matchs')==0,     options.matchs=[]; end
        if isfield(options,'CF_match')==0,     options.CF_match='D'; end
end