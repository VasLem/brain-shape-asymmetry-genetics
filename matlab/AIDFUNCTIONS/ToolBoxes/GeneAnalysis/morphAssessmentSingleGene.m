function ass = morphAssessmentSingleGene(IndVar,shape,RefScan,BL,name)
         res = morphingSingleGenes(IndVar,shape,RefScan,BL);
         ass = assessment;
         ass.Scan = clone(res(1).scan);
         ass.Scan.Tag = [name '_M'];
         ass.Norm = clone(res(2).scan);
         ass.Norm.Tag = [name '_P'];
         ass.DistanceRange = [0 0.1];
         ass.Significance = 6;
         update(ass);
         ass.Tag = name;
end