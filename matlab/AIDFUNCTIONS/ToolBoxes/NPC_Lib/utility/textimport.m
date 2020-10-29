function [D,data,names]=textimport(varargin);
% [D,data,names]=textimport(FILE,DELIMITER,HEADER)
%     FILE: string defining the file to read from. Default directory is pwd.
%           Default extension is 'txt'.
%     DELIMITER: in string format. '\t' by default
%     HEADER: if HEADER is present: 1 (by default), otherwise: 0
%
%see also csvexport, xlsimport
%Livio Finos
%e-mail: livio@stat.unipd.it


if length(varargin)<2
    varargin{2}='\t';
end
if length(varargin)<3
    varargin{3}=1;
end


dlm=varargin{2};
filename=varargin{1};
if varargin{3}==1
    [data]=dlmread(filename,dlm,1,0);
    [n m]=size(data);
    [names]=textread(filename,'%s',m,'delimiter',dlm); %,'emptyvalue','NONAME'
else
    [data]=dlmread(filename,dlm,0,0);
    [n m]=size(data);
    for i=1:m
        names{i}=['V' num2str(i)];
    end
end




for i=1:m
    D(i).name=names(i);
    D(i).vals=data(:,i);
    D(i).code=cell(0,3);
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
