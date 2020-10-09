function out = createG2FModel(input,GENESNAMES,RIPGENES,GENES,SEX,RIPS,ANC,RIPA,Model)
         % assembling required sex info
         out.Findex = find(SEX==1);
         out.Mindex = find(SEX==-1);
         out.RIPS = RIPS;
         % assembling required ancestry info
         out.RIPA = RIPA;
         [out.AncP,out.AncS] = polyfit(ANC,RIPA,1);
         % assembling required gene info
         out.Genes = input.GENENAMES;
         out.nrGenes = length(out.Genes);
         out.nrSamples = length(SEX);
         out.RIPG = zeros(out.nrSamples,out.nrGenes);
         out.GENES = zeros(out.nrSamples,out.nrGenes);
         out.HETERO = zeros(1,out.nrGenes);
         out.ANOVAheadings = input.ANOVAHeadings;
         out.ANOVA = input.ANOVA;
         out.Model = clone(Model);
         out.Shape = reconstructTraining(out.Model)';
         for i=1:1:out.nrGenes
             index = find(strcmp(out.Genes{i},GENESNAMES));
             out.RIPG(:,i) = RIPGENES(:,index);
             out.GENES(:,i) = GENES(:,index);
             FAA = out.ANOVA(i,5);
             FBB = out.ANOVA(i,7);
             if FAA<=FBB
                out.HETERO(i) = -1;
             else
                out.HETERO(i) = 1;
             end          
         end
end