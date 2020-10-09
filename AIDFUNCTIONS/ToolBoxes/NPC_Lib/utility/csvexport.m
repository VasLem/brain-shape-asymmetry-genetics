function []=csvexport(varargin);
% []=csvexport(FILE,D)
% write the data D in to the worksheet SHEET (by default 'Data'), in the Excel
% file specified in FILE (by default 'NPCData.xls').
% The data in the xls sheet are a nxm matrix (n number of observations and
% m variables) with variables name in the first row (does not use accented 
% characters). see also example_xlsimport.xls in the DIR\NPClib\examples folder. 
%
%     INPUT PARAMETERS (See xlswrite.m for other options):
%     FILE: string defining the file to read from. Default directory is pwd.
%           Default extension is 'xls'.
%     D:    a struct array containing the data in NPC format.
%
%
%Livio Finos
%e-mail: livio@stat.unipd.it


D=varargin{2};

fid = fopen(varargin{1},'wb');
C=[];
for i=1:length(D)
    C=[C char(D(i).name) ','];
end
C=[C '\n'];

for j=1:length(D(1).vals)
    
for i=1:length(D)
    if isnan(D(i).vals(j))
        C=[C  ','];
    else
        if size(D(i).code,2)
            lb=find([D(i).code{:,2}]==D(i).vals(j));
            if not(isempty(lb))
                C=[C D(i).code{lb(1),1} ','];
            else
                C=[C num2str(D(i).vals(j)) ','];
            end
        else
            C=[C num2str(D(i).vals(j)) ','];
        end
    end
end

C=[C '\n'];
end

fprintf(fid,C);
fclose(fid);
fprintf('\n\tData have been exported on file: %10s \n', varargin{1});
    fprintf('\n\tSummary dataset:\n%11.0f\t Rows\n%11.0f\t Columns\n\n',size([D.vals]))