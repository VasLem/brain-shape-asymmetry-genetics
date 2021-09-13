function [reducedTemplate, reducedLH,reducedRH, landmarksIndices] = preprocessSymmetry(template, LH, RH, reduction_rate, subsampling_rate)
%PREPROCESSHEMISPHERES Preprocess LH and RH for symmetry analysis
%1. Mirrors RH on LH by x axis reversion.
%2. Reduces landmarks number provided LH and RH if reduction_rate is given using getDownsampledLandmarksIndices function.
%3. Aligns them using Procrustes' Analysis GeneralizedProcrustesAnalysis function
%4. Subsamples them using subsampling_rate, if it has been provided
%Template is the shape3D object corresponding to LH
disp("Mirroring RH to LH..");
RH(:,1,:,:) = -1*RH(:,1,:,:);
if (nargin > 3) && (reduction_rate<1)
    disp("Retrieving Indices for Downsampling MRI Image..")
    [landmarksIndices, reducedFaces, ~]  = getDownsampledLandmarksIndices(template,reduction_rate,true);
else
    landmarksIndices = 1:size(LH,1);
    reducedFaces = template.Faces;
end
reducedLH = LH(landmarksIndices,:, :);
reducedRH = RH(landmarksIndices,:, :);
reducedTemplate = shape3D;
reducedTemplate.Vertices = template.Vertices(landmarksIndices, :) ;
reducedTemplate.Faces = reducedFaces;

disp("Applying Procrustes Analysis..")
[AlignedShapes,~,~] = GeneralizedProcrustesAnalysis(cat(3, reducedLH, reducedRH), reducedTemplate,3,true,false,true,false);
reducedLH = AlignedShapes(:,:,1:size(reducedLH,3));
reducedRH = AlignedShapes(:,:,size(reducedLH,3)+1:end);

len = size(reducedLH,3);
if (nargin > 4) && (subsampling_rate<1)
    disp("Applying subsampling")
    step = 1 / subsampling_rate;
    samplesIndices = 1:step:len;
else
    samplesIndices = 1:len;
end
reducedLH = reducedLH(:, :, samplesIndices);
reducedRH = reducedRH(:, :, samplesIndices);
end

