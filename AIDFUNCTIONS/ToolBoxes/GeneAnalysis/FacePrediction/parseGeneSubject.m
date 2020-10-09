function out = parseGeneSubject(input)
         [~,~,raw] = xlsread(input);
         out.Sex = raw{1,2};
         out.Anc = cell2mat(raw(2,2));
         out.Genes = raw(3:end,1);
         out.GenoTypes = cell2mat(raw(3:end,2));
         out.nrGenes = length(out.Genes);
end