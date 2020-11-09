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

da_values <- rep(sym$DA.component,each=Ns)
fa_values <- rep(sym$FA.component,each=Ns)
template <- GPA_DATA_INPUT$strTemplate
