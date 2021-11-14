function [out, score] = computeAmmiModel(Shapes,  numComponents)
arguments
Shapes
numComponents=0
end
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
flucMat = sideDif - mean(sideDif,3);
flucMat = reshape( permute(flucMat, [2,1,3]), [],size(flucMat,3));
if numComponents==0
disp("..Selecting optimal number of Principal Components to retain (Parallel Analysis)..")
[latent, ~, latentHigh] = parallelAnalysis(flucMat, 100,0.05);
numComponents = sum(latent > latentHigh');
end
disp(['..Selected number of components: ' num2str(numComponents)]);
disp("..Applying PCA reconstruction..")
[coeff, score, ~] =  pca(flucMat,'Centered',false, 'NumComponents', numComponents);
rec = coeff * score';
rec = mean(rec,1);
rec = sqrt(sum(reshape(rec,3,[]).^2, 1)');
score = squeeze(sqrt(sum(reshape(score, 3, [], numComponents).^2, 1)));
F = rec;
out.LM.I = I';
out.LM.D = D';
out.LM.F = F';
end
