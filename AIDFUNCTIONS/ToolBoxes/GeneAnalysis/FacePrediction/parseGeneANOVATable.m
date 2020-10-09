function out = parseGeneANOVATable(input)
         [~,~,raw] = xlsread(input);
         out.GENENAMES = raw(2:end,1);
         out.ANOVAHeadings = raw(1,2:end);
         out.ANOVA = cell2mat(raw(2:end,2:end));
end