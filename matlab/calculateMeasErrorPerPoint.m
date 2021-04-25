addpath(genpath('AIDFUNCTIONS'));
addpath(genpath('fieldtrip'));
DATA_DIR = '../SAMPLE_DATA/';
DATA = cell(1,2);
Render = cell(1,2);
phenopath = [DATA_DIR, 'IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/'];
Regions = {'LH' 'RH'};
r=1;
regphenopath = [phenopath Regions{r} '/'];
%r=2
brainSurface = load([regphenopath 'RENDERMATERIAL.mat']);
scale = imread('../SAMPLE_DATA/MeasError/scale0to4pc.png');
scale_colors = uint8(squeeze(mean(scale,1)));
scale_values = 1:(3/size(scale_colors,1)):4;
scale_values = scale_values(2:end);

%%
front = imread('../SAMPLE_DATA/MeasError/front.png');
back = imread('../SAMPLE_DATA/MeasError/back.png');
%%

scan = brainSurface.RefScan;
map = combvec(1:5:250,1:5:250,1:5:250);
map = map(:,1:scan.nFaces)'/256;
% ,scan.nFaces,2,'Replace',false)'/256;
% r.VertexRGB = map(:,:);
scan.ColorMode = "texture";
viewer = scan.viewer();
scan.PatchHandle.FaceColor ='flat';
scan.PatchHandle.FaceVertexCData = map;

%%
view(gca(),brainSurface.viewval(1),0);
exportgraphics(gcf(),'labeled_front.png', 'resolution',300);
view(gca(),-brainSurface.viewval(1),0);
exportgraphics(gcf(),'labeled_back.png', 'resolution',300);
%%
close all
ref = imread('labeled_front.png');
imshow(ref)

frontVerticesValues = mapToRef(front, ref,scale_colors,scale_values,map,scan);
%%
ref = imread('labeled_back.png');
imshow(ref)
backVerticesValues = mapToRef(back, ref,scale_colors,scale_values,map,scan);
%%
outputValues = frontVerticesValues;
outputValues(outputValues==0) = backVerticesValues(outputValues==0) ;
%%
adjacency = scan.Adjacency;

outputValues = reconstruct(outputValues,adjacency);

%%


%%

outSurface = load([regphenopath 'RENDERMATERIAL.mat']);
scan = outSurface.RefScan;
scan.VertexValue    = outputValues;

scan.ColorMode = "Indexed";
scan.Material = 'Dull';
scan.ViewMode = 'solid';
v =scan.viewer();

colormap(scale_colors);

%%
view(gca(),outSurface.viewval(1),0);
%%
view(gca(),-outSurface.viewval(1),0);
%%
percent_difference_test_retest = outputValues;
save('../SAMPLE_DATA/MeasError/percent_difference_test_retest','percent_difference_test_retest');

%%
function reconstructed = reconstruct(values,adjacency)
values(values==0) = NaN;
reconstructed = values;
pdim = size(values,1);
for l = 1:50
    for c= 1:pdim
        if isnan(reconstructed(c))
           f = reconstructed((adjacency(:,c)>0) & ~ squeeze(isnan(reconstructed)));
           if ~isempty(f)
               reconstructed(c) =  mean(f,1);
           end
        end
    end
    if sum(isnan(reconstructed))==0
        break
    end
end
end

function verticesValues=mapToRef(to_cmp, ref,scale_colors,scale_values, map,scan)
bw_ref = (~(ref(:,:,1) == 255 & ref(:,:,2) == 255 & ref(:,:,3) ==255));
ref = ref .* uint8(bw_ref) ;
bw = (~(to_cmp(:,:,1) <= 10 & to_cmp(:,:,2) <=10 & to_cmp(:,:,3) <=10));
l = bwlabel(bw);
counts=hist(l(l~=0),1:10);
bw = l==argmax(counts);
imshow(bw)
%%
moving = imresize(double(bw), size(bw_ref), 'nearest')>0;
fixed = bw_ref;
%%
[D,movingReg] = imregdemons(moving,fixed,[500 400 200],...
    'AccumulatedFieldSmoothing',1.3);
%%
imshowpair(bw_ref,movingReg,'montage');
%%
ret = imwarp( imresize(to_cmp, size(bw_ref), 'nearest'),D) .* uint8(moving);
imshowpair(ref .* uint8(bw_ref) ,ret,'montage')


%%
% Keeping only colors existing in scale
m = zeros(size(ret,[1,2]));
ret = double(ret);
scale_colors = double(scale_colors);
t = 10;
for i=1:size(scale_colors, 1)
    m = m + (( (abs(ret(:,:,1) - scale_colors(i,1)) < t)...
        & (abs(ret(:,:,2) -scale_colors(i,2))<t)...
        & (abs(ret(:,:,3) - scale_colors(i,3))<t)) );
end
ret = uint8(ret);
scale_colors = uint8(scale_colors);
filtered = uint8(m>0).* uint8(ret);
imshow(filtered);
%%

maskToUse = squeeze(sum(filtered > 0,3)>0);
refToUse = ref .* uint8(maskToUse);
%%
% Mapping to Annotation Map
t = 4;
mappedColors = zeros(size(map));
u8map = reshape(double(uint8(256 * map)),[1,size(map,1),3]);
refToUse = double(refToUse);
parfor i=1:size(map,1)
    m = max(abs(refToUse - u8map(:,i,:)),[],3)<t;
    norm = sum(m,'all');
    if norm<20
        continue
    end
    f = filtered .* uint8(m);
    f = squeeze(sum(f,[1,2]));
    f = uint8(double(f)/norm);
    mappedColors(i,:) = f;
end
%%
% Mapping to scale
t = 10;
scale_colors = double(scale_colors);
mappedColors = double(mappedColors);
mappedValues = zeros(size(mappedColors,1),1);
for i=1:size(mappedColors,1)
    if all(mappedColors(i,:)==0)
        continue
    end
    for j=1:size(scale_colors,1)
        m = max(abs(mappedColors(i,:) - scale_colors(j,:)),[],2) < t;
        if m
            mappedValues(i) = scale_values(j);
            break
        end
    end
end
%%
% Mapping Faces values to Vertices
verticesValues = zeros(scan.nVertices,1);
for i=1:size(mappedValues,1)
    if mappedValues(i) == 0
        continue
    end
    verticesValues(scan.Faces(i,1)) = mappedValues(i);
    verticesValues(scan.Faces(i,2)) = mappedValues(i);
    verticesValues(scan.Faces(i,3)) = mappedValues(i);
end
%%
end


