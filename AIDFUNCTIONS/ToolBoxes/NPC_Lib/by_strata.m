function [P T options]=by_strata(strata, namef, Y,varargin);

if isstruct(strata)
    name=[strata.name];
    Dstrata=strata;
    strata=[strata.vals];
elseif ischar(strata)
    name=strata;
    [strata, Dstrata]=v({strata});
elseif iscell(strata)
    name=strata{1};
    [strata, Dstrata]=v(strata);
else
    strata=strata(:);
    name='strata';
    Dstrata.code=cell(0,3);
end

if isempty(strata)
    strataname=cell(1,max([Dstrata.code{:,2}]));
    strataname([Dstrata.code{:,2}])=[Dstrata.code(:,1)];
else
    strataname=cell(1,max(strata));
    levs=unique(strata)';
    for i=1:length(levs), strataname{i}=num2str(levs(i)); end
end

if ischar(Y) | iscell(Y)
    [no Y]=v(Y);
end

if not(strcmpi(namef,'NP_1s'))
    X=varargin{1};
    varargin=varargin(2:end);
    if isnumeric(X)& isvector(X)
        X=X(:);    
    elseif ischar(X) | iscell(X)
        [no X]=v(X);
    end
    if strcmpi(namef,'NP_rho')    
        Z=varargin{1};
        varargin=varargin(2:end);
        if isnumeric(Z)& isvector(Z)
            Z=Z(:);    
        elseif ischar(Z) | iscell(Z)
            [no Z]=v(Z);
        end
    end
end


levs=unique(strata);
%if strataname num2str(levs(i))

for i=1:length(levs)
    fprintf(['\nStrata=' strataname{i}]    )
    [y]=getvalues(Y,strata==levs(i));
    if strcmpi(namef,'NP_1s')
        eval(['[P(i,:,:,:) T(i,:,:,:) options]=' namef '(y, varargin{:} );'])
    else
        [x]=getvalues(X,strata==levs(i));
        if strcmpi(namef,'NP_rho')    
            [z]=getvalues(Z,strata==levs(i));
            eval(['[P(i,:,:,:,:) T(i,:,:,:,:) options]=' namef '(y,x,z,varargin{:} );'])
        else
            eval(['[P(i,:,:,:) T(i,:,:,:) options]=' namef '(y,x,varargin{:} );'])
        end
    end
end
p=P(:,end,:,:); p=shiftdim(p,1);
P=shiftdim(P,1);
T=shiftdim(T,1);
    
ndims=length(size(P));

options.p.raw=p;
options.labels.dimslabel{ndims}={name};
options.labels.dims{ndims}=strataname;

function [w]=getvalues(W, id)
if isstruct(W)
    w=W;
    for j=1:length(W)
        w(j).vals=W(j).vals(id);
    end
else
    w=W(id,:,:);
end
