function out = computeAmmiModel(Shapes)
inputShapes = Shapes;
n = size(inputShapes,3)/2;
LH = inputShapes(:,:,1:n);
RH =  inputShapes(:,:,n+1:2*n);
avgShapes = (LH + RH)/2;
sideDif =LH  - RH ;
%%
disp("Computing Directional Asymmetry..")
D = sqrt(sum(mean(sideDif,3).^2,2));
avgShape = mean(avgShapes,3);
disp("Computing Individual Asymmetry..")
indDif = avgShapes - avgShape ;
I = mean(sqrt(sum(indDif.^2,2)),3);
%%
disp("Computing Fluctuating Asymmetry..")
flucMat = LH - RH + avgShape;
flucMat= squeeze(sqrt(sum((flucMat).^2,2)));
disp("..Selecting optimal number of Principal Components to retain (Parallel Analysis)..")
[latent, ~, latentHigh]  = parallelAnalysis(flucMat', 100,0.05);
numComponents = sum(latent > latentHigh');
disp(["..Selected number of components: " num2str(numComponents)]);
disp("..Applying PCA reconstruction..")
[coeff, score, ~] =  pca(flucMat','Centered',false);
rec = coeff(:, 1:numComponents) * score(:, 1:numComponents)' ;
 F = mean(rec,2);
out.LM.I = I';
out.LM.D = D';
out.LM.F = F';
out.Raw.F = rec;
end

