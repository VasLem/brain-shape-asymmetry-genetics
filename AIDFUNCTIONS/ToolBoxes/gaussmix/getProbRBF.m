function ProbRBF = getProbRBF(UV, prob)
         % UV = newCUV;
         pl.Location = UV;
         pl.Value = prob';
         %pl = fastrbf_unique(pl,'messages',0);
         pl = uniquePL(pl);
         ProbRBF = fastrbf_fit(pl,0.01,'reduce','messages',0);
end