function [out,model,v] = evalGroupModel(DATA,Template,display)
         if nargin<3, display = false; end
         N = length(DATA);
         Shapes = zeros(Template.nVertices,3,N);
         for i=1:1:N
            Shapes(:,:,i) = DATA{i}.MappedShape.Vertices;
         end
         % align all shapes
         [AlignedShapes,AverageShape,~,v] = GeneralizedProcrustesAnalysis(Shapes,Template,3,true,false,display);
         % center around the average shape
         CenteredShapes = AlignedShapes-repmat(AverageShape.Vertices,1,1,N); 
         % redimension shapes
         CenteredShapes = reshape(CenteredShapes,Template.nVertices*3,N);
         warning off;
         [COEFF, SCORE, LATENT, ~, EXPLAINED] = pca(CenteredShapes','Economy',true,'Algorithm','svd','Centered',false);
         warning on;
         model.EigVec = COEFF;
         model.EigVal = LATENT;
         model.Scores = SCORE;
         model.VarExp = EXPLAINED;        
         out.top10 = sum(EXPLAINED(1:10));
         out.top5 = sum(EXPLAINED(1:5));
         count = 0;total = 0;
         while total<95
             count = count+1;
             total = total+EXPLAINED(count);
         end
         out.Var98 = count;
         out.CumVar = cumsum(EXPLAINED);
         if display, figure;plot(out.CumVar);end
end