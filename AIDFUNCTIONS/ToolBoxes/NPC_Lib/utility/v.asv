function [vals, Dreduct]=v(varargin);
%vals=v(idcolum)
%vals=v('varname1','varname2',...); and/or
%vals=v(columnnumbers,'varname1',...);
%after use of function reminD.m
%
% 
%load germina , reminD(Dgermina)
%stesso risultato facendo: [Dgermina(2:4).vals] oppure v(2:4)
% oppure v({Germinated','Weight','Surface'}) oppure
% v(Germinated','Weight','Surface')
%
%possono essere selezionate più variabili con un solo argomento:
% es: v('T') richiama tutte le variabili che cominciano per 'T'
%
%Livio Finos
%e-mail: livio@stat.unipd.it

global NPCD__
vals=[];
Dreduct=NPCD__(1:-1);

if isempty(NPCD__)
    fprintf('\n No data are reminDed by NPClib.\n Make use of function reminD.m before v.m')
else
    
    
    temp=cell(0,0);
for i=1:length(varargin) 
    if iscell(varargin{i})
        for j=1:length(varargin{i})
            temp{end+1}=varargin{i}{j};
        end
    else
        temp{end+1}=varargin{i};
    end
end
varargin=temp;

Dreduct=NPCD__([]);

for i=1:length(varargin)
    if isnumeric(varargin{i})
        vals=[vals NPCD__(varargin{i}).vals];
        Dreduct(end+[1:length(varargin{i})])=NPCD__(varargin{i});
    else
        quo=strmatch(varargin{i},[NPCD__.name],'exact');
       if isempty(quo) quo=strmatch(varargin{i},[NPCD__.name]); end
       
        vals=[vals NPCD__(quo).vals];
        if not(isempty(quo))
            Dreduct(end+(1:length(quo)))=NPCD__(quo);
        end
        
    end
end
end
