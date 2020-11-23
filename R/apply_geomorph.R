library(knitr)
library(rgl)
knit_hooks$set(webgl = hook_webgl)
require(geomorph)
library(R.matlab)
# install.packages("tidyverse")
library(tidyr)
library(scatterplot3d) # load

GPA_DATA <- readMat("../SAMPLE_DATA/m2r.mat")
nSamples <- GPA_DATA$nSamples
Ns <- GPA_DATA$Ns
nReps <- GPA_DATA$nRep


# center
# centers a matrix faster than scale()
# used in various functions where mean-centering is required
center <- function(x){
  if(is.vector(x)) x - mean(x) else {
    x <- as.matrix(x)
    dims <- dim(x)
    fast.center(x, dims[1], dims[2])
  }
}

fast.center <- function(x, n, p){
  m <- colMeans(x)
  x - rep.int(m, rep_len(n, p))
}

fast.scale <- function(x, n, p){
  if(p > 1) {
    x <- fast.center(x, n, p)
    scale <- apply(x, 2, sd)
    x / rep.int(scale, rep_len(n, p))
  } else {
    x <- x - mean(x)
    x/sd(x)
  }
}

# csize
# calculates centroid size
# digitsurface
csize <- function(x) sqrt(sum(center(as.matrix(x))^2))

LH <- GPA_DATA$LH
RH <- GPA_DATA$RH
template <- GPA_DATA$strTemplate[[2]]
dim(template)

allShapes <- array(c(template, LH, RH), dim = c(dim(LH)[1], dim(LH)[2], 1 + dim(LH)[3] + dim(RH)[3]))

Y.gpa = gpagen(allShapes, Proj=FALSE, PrinAxes=FALSE,max.iter=3)
templateAligned = Y.gpa$coords[,,1];
LHAligned=Y.gpa$coords[,,2:(dim(LH)[3]+1)]
RHAligned=Y.gpa$coords[,,-(1:(dim(LH)[3] + 1))]

# scalingFactor = csize(template)/csize(TemplateAligned);
# scaledTemplateAligned = TemplateAligned * scalingFactor;

# translationFactor = colMeans(sampledTemplate)- colMeans(scaledTemplateAligned); 
# recSampledTemplate = scaledTemplateAligned + translationFactor;


writeMat("../SAMPLE_DATA/r2m.mat",LHAligned=LHAligned,RHAligned=RHAligned, TemplateAligned=templateAligned)


ind <- c(0, rep(1:nSamples, 2*nReps))
rep <- c(1 ,rep(rep(1:3,each=nSamples),2))
side <- c(1, c(rep(1,nSamples*nReps), rep(2,nSamples*nReps)))

ind <- c(0, rep(1:nSamples, 2))
rep <- c(1 ,rep(rep(1,each=nSamples),2))
side <- c(1, c(rep(1,nSamples*1), rep(2,nSamples*1)))



sym <- bilat.symmetry(Y.gpa, ind=ind,side=side,rep=rep)

plot(sym, warpgrids = TRUE)
dim(sym$DA.component)
da_values <- rowSums((sym$DA.component[,,1]-sym$DA.component[,,2])^2)
da_values <- rep(da_values,each=Ns)
fa_values <- rowSums((sym$FA.component[,,1]-sym$FA.component[,,2])^2)
fa_values <- rep(fa_values,each=Ns)
template <- GPA_DATA_INPUT$strTemplate

vertices <- rbind(t(as.matrix(template[[2]])),1)
indices <- template[[3]]



# center(template)
# templateCenter <- mean(template, 1)
# templateCenter
# 
# 
# 
# 
# material <- fa_values[1:dim(vertices)[2]]
# 
# 
# 
# material <- as.integer(255 * (material-min(material))/(max(material)-min(material)))
# paletteFunc <- colorRampPalette(c("blue", "red"));
# 
# palette <- paletteFunc(255);
# 
# mesh1 <- tmesh3d(vertices = vertices, indices = indices, homogeneous = TRUE, 
#                  material = list(palette[material])
#                  )
# 
# shade3d(mesh1, meshColor = "vertices")
# 
# 
# wire3d(mesh1)
