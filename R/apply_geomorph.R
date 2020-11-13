library(knitr)
library(rgl)
knit_hooks$set(webgl = hook_webgl)
require(geomorph)
library(R.matlab)
library(abind)

GPA_DATA_INPUT <- readMat("../SAMPLE_DATA/input_to_r.mat")


LH <- GPA_DATA_INPUT$LHReps
RH <- GPA_DATA_INPUT$RHReps
allShapes <- abind(LH, RH,along=3)
nSamples <- GPA_DATA_INPUT$nSamples
Ns <- GPA_DATA_INPUT$Ns
nReps <- GPA_DATA_INPUT$nRep


ind <- rep(1:nSamples, 2*nReps)
rep <- rep(rep(1:3,each=nSamples),2)
side <- c(rep(1,nSamples*nReps), rep(2,nSamples*nReps))

sym <- bilat.symmetry(allShapes, ind=ind,side=side,rep=rep)

plot(sym, warpgrids = TRUE)
dim(sym$DA.component)
da_values <- rowSums((sym$DA.component[,,1]-sym$DA.component[,,2])^2)
da_values <- rep(da_values,each=Ns)
fa_values <- rowSums((sym$FA.component[,,1]-sym$FA.component[,,2])^2)
fa_values <- rep(fa_values,each=Ns)
template <- GPA_DATA_INPUT$strTemplate

vertices <- rbind(t(as.matrix(template[[2]])),1)
indices <- template[[3]]


material <- fa_values[1:dim(vertices)[2]]



material <- as.integer(255 * (material-min(material))/(max(material)-min(material)))
paletteFunc <- colorRampPalette(c("blue", "red"));

palette <- paletteFunc(255);

mesh1 <- tmesh3d(vertices = vertices, indices = indices, homogeneous = TRUE, 
                 material = list(palette[material])
                 )

shade3d(mesh1, meshColor = "vertices")


wire3d(mesh1)
