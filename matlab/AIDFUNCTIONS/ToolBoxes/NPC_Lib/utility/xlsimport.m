function [D,data,code]=xlsimport(varargin);
% [D,data,code]=xlsimport(FILE,SHEET,RANGE,MODE) 
% reads the data specified in RANGE from the worksheet SHEET, in the Excel
% file specified in FILE. 
% The data in the xls sheet are a nxm matrix (n number of observations and
% m variables) with variables name in the first row (does not use accented 
% characters). see also example_xlsimport.xls in the DIR\NPClib\examples folder. 
%
%     INPUT PARAMETERS (See xlsread.m for other options):
%     FILE: string defining the file to read from. Default directory is pwd.
%           Default extension is 'xls'. See NOTE 1.
%     SHEET: string defining worksheet name in workbook FILE.
%            double scalar defining worksheet index in workbook FILE.
%     RANGE: string defining the data range in a worksheet. See NOTE 2.
%     MODE: string enforcing basic import mode. Valid value = 'basic'.
%
%     OUTPUT PARAMETERS:
%     D:    a struct array containing the data and ready to be used in
%     NPClib.
%     D.vals
%     D.name
%     D.code
%     data: the nXm matrix of numerical the data.
%     code: the structure array containing recoding information for string obervations.
%
%
%Livio Finos
%e-mail: livio@stat.unipd.it




[data labels]=xlsread(varargin{:});

text=labels(2:end,:);
labels=labels(1,:);


check=1;
if size(text,1)>0
while check>0
    if length(strmatch('',text(:,check),'exact'))==0
        data=[NaN.*ones(size(data,1),1) data];
        check=check+1;
        if check>size(text,2)
            check=0;
        end
    else
        check=0;
    end
end
end

if size(data,2)<length(labels)
    data(:,end+1:length(labels))=NaN;
end

if size(data,1)>0
if (size(text,1)>0) & (length(find(isnan(data(1,:))))==length(data(1,:)))
    data=data(2:end,:);
end
end
    

for i=1:length(labels)
    code(i).labels=labels(1,i);
    code(i).etich=unique(text(:,i));
    code(i).etich=code(i).etich(setdiff(1:length(code(i).etich),strmatch('NaN',code(i).etich,'exact')));
    code(i).etich=code(i).etich(setdiff(1:length(code(i).etich),strmatch('',code(i).etich,'exact')));
    temp2=unique(data(:,i));
    temp2=temp2(find(not(isnan(temp2))));
    temp=setdiff(1:1000,temp2);    
    for j=1:length(code(i).etich)
        quali=strmatch(code(i).etich(j),text(:,i),'exact');
        data(quali,i)=temp(j);
        code(i).etich{j,2}=temp(j);
        code(i).etich{j,3}=labels{i};
    end 
    j=size(code(i).etich,1);
    
    if sum(isnan(temp2))<length(temp2) & not(isempty(code(i).etich))
        for j2=1:length(temp2)
            code(i).etich{j+j2,1}=num2str(temp2(j2));
            code(i).etich{j+j2,2}=temp2(j2);
            code(i).etich{j+j2,3}=labels{i};
        end
    end
    if not(isempty(code(i).etich))
    for j=1:size(code(i).etich,1)
        quali=strmatch(code(i).etich(j),text(:,i),'exact');
        data(quali,i)=temp(j);
        code(i).etich{j,2}=temp(j);
        code(i).etich{j,3}=labels{i};
    end    
    end
end

for i=1:length(labels)
   D(i).name=code(i).labels;
   D(i).code=code(i).etich;
   D(i).vals=data(:,i);
end

if 1 %not(OUT==0)
    fprintf('\n\tSummary dataset:\n%11.0f\t Observations\n%11.0f\t Variables\n\t________________________________________',size([D.vals]))
    type={'cathegorical','quantitative'};
    fprintf('\n %20s\t%20s','Name', 'Type')
    for i=1:length(D)
        fprintf('\n %20s\t%20s',D(i).name{1}(1:min(length(D(i).name{1}),20)), type{isempty(D(i).code)+1 })
    end
    fprintf('\n\t________________________________________')
    for i=1:length(D)
        if not(isempty(D(i).code))
            fprintf('\n\t Codes of variable %s',D(i).name{1})
            fprintf('\n %10s\t %10s','Label', 'Value')
            for j=1:size(D(i).code(:,1:2),1)
                fprintf('\n %10s\t %10s',D(i).code{j,1},num2str(D(i).code{j,2}) )
            end
            fprintf('\n')
        end
    end
    fprintf('\n')
end
