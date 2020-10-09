function score = myLocalSolverFun(obj,X)
         score = getScore(obj.DemiModelSet,X')...
                 + obj.FacePenalty*sqrt(sum(((X-obj.FaceModel.AvgCoeff')./obj.FaceModel.EigStd').^2)); 
end