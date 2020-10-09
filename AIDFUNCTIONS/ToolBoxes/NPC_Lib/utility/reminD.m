function []=reminD(D,OUT);
% reminD(D,OUT); 
% The dataset to remind is a dataset in NPClib format, as imported by
% xlsimport or textimport.
% if OUT=0 descriptives are not printed.
%
% see also: v, xlsimport, textimport.
%
% NPClib 
%Livio Finos
%e-mail: livio@stat.unipd.it

global NPCD__
NPCD__=D;

if nargin==1
    OUT=1;
end
if not(OUT==0)
    fprintf('\n\t___________________________\n\tReminD dataset:\n%11.0f\t Observations\n%11.0f\t Variables\n\t___________________________\n',size([D.vals]))
end
