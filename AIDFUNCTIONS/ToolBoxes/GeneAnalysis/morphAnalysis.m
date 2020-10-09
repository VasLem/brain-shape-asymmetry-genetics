function out = morphAnalysis(A,shape,RefScan,scale,CondValues)
         % MorphCreation
           %out.morphs = createMorphs(A,shape,RefScan,scale,CondValues);
           [out.morphs,out.STATS] = createMorphsReducedModel(A,shape,RefScan,scale,CondValues);
         % Compute image MorphAnalysis  
           out.comparison = compareMorphs(out.morphs(1).scan,out.morphs(2).scan);
         % compute effect and effect-size
         
end