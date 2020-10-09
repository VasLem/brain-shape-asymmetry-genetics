function [GV,codings,codingsV] = transCodeGene(G)
%           G = cell(1,5);
%           G{1} = 'AA';
%           G{2} = 'AB';
%           G{3} = 'BA';
%           G{4} = 'BB'; 
          nrSamples = length(G);
          GV = nan*ones(1,nrSamples);
          % Detecting good ones
          index = [];
          badindex = [];
          for i=1:1:nrSamples
              if ~isempty(G{i})&~isnan(G{i})
                 index = [index i]; %#ok<AGROW>
              else
                 badindex = [badindex i];
              end
          end
          Gred = G(index);
          codings = unique(Gred);
          nrCodes = length(codings);
          % Assingings codingsV
          codingsV = zeros(1,nrCodes);
          for i=1:1:nrCodes
              tmp = codings{i};
              if strcmp(tmp(1),tmp(2))
                  codingsV(i) = 1;
              end
          end
          ind = find(codingsV);
          codingsV(ind(1)) = -1;
          % assinging codes
          for i=1:1:length(index)
              sample = G{index(i)};
              tmp = strcmp(codings,sample);
              ind = tmp;
              GV(index(i)) = codingsV(ind);
          end      
end