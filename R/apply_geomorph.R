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
LH <- GPA_DATA$LH
RH <- GPA_DATA$RH
template <- GPA_DATA$strReducedTemplate[[2]]
templateOriginal <- GPA_DATA$strTemplate[[2]]

# center
# centers a matrix faster than scale()
# used in various functions where mean-centering is required
center <- function(x) {
  if (is.vector(x))
    x - mean(x)
  else {
    x <- as.matrix(x)
    dims <- dim(x)
    fast.center(x, dims[1], dims[2])
  }
}

fast.center <- function(x, n, p) {
  m <- colMeans(x)
  x - rep.int(m, rep_len(n, p))
}

fast.scale <- function(x, n, p) {
  if (p > 1) {
    x <- fast.center(x, n, p)
    scale <- apply(x, 2, sd)
    x / rep.int(scale, rep_len(n, p))
  } else {
    x <- x - mean(x)
    x / sd(x)
  }
}

# csize
# calculates centroid size
# digitsurface
csize <- function(x)
  sqrt(sum(center(as.matrix(x)) ^ 2))


dim(template)
library(abind)


allShapes <- abind(template, LH, RH, along = 3)
allShapes.gpa = gpagen(allShapes,
                       Proj = FALSE,
                       PrinAxes = FALSE,
                       max.iter = 3)
allShapes.gpa$template = allShapes.gpa$coords[, , 1]

allShapes.gpa$LH = allShapes.gpa$coords[, , 2:(dim(LH)[3] + 1)]
allShapes.gpa$RH = allShapes.gpa$coords[, , -(1:(dim(LH)[3] + 1))]

allShapes.gpa$LHReps = allShapes.gpa$LH

function(x) {
  m <- mean(x)
  rnorm(some_number, m, 0.1*m)
}

mat <-allShapes.gpa$LH
createNoisySample<- function(mat){
  dim(mat) <- c(298 * 3 , 100)
  w <-apply(mat,2,var)
  m <- matrix(rnorm(prod(dim(mat))), dim(mat)) * w
  noisy <- mat + m * 0.2
  dim(noisy) <- c(298,3,100)
  noisy
}

createNoisyHemispheres <- function(hemi){
  ret <- c()
  for (rep in 1:nReps){
    ret <- abind(ret, createNoisySample(hemi), along=3);
  }
  ret
}
allShapes.gpa$noisyAllShapes = abind(allShapes.gpa$template,
  createNoisyHemispheres(allShapes.gpa$LH), createNoisyHemispheres(allShapes.gpa$RH),along=3)
dim(allShapes.gpa$noisyAllShapes)


writeMat(
  "../SAMPLE_DATA/r2m.mat",
  LHAligned = allShapes.gpa$LH,
  RHAligned = allShapes.gpa$RH,
  TemplateAligned = allShapes.gpa$template
)


# ind <- c(0, rep(1:nSamples, 2*nReps))
# rep <- c(1 ,rep(rep(1:3,each=nSamples),2))
# side <- c(1, c(rep(1,nSamples*nReps), rep(2,nSamples*nReps)))

ind <- c(0, rep(1:nSamples, 2*nReps))
rep <- c(1 , rep(rep(1:nReps, each = 2 * nSamples)))
side <- c(1, rep(c(rep(1, nSamples * 1), rep(2, nSamples * 1)),nReps))

nPoints <- dim(allShapes.gpa$noisyAllShapes)[1]
toPlot <- seq(1, nPoints, 10)
nToPlot <- length(toPlot)

sym <- bilat.symmetry(allShapes.gpa$noisyAllShapes,
                      ind = ind,
                      side = side,
                      rep=rep,
                      iter = nSamples)
summary(sym)
plot(sym, warpgrids = TRUE)
dim(sym$DA.component)
da_values <-
  rowSums((sym$DA.component[, , 1] - sym$DA.component[, , 2]) ^ 2)
# da_values <- rep(da_values, each = Ns)

fa_values <-
  rowSums((sym$FA.component[, , 1] - sym$FA.component[, , 2]) ^ 2)
# fa_values <- rep(fa_values, each = Ns)

ia_values <-
  rowSums((sym$IA.component[, , 1] - sym$IA.component[, , 2]) ^ 2)
# fa_values <- rep(fa_values, each = Ns)  


template <- GPA_DATA$strTemplate

representComponent <- function(component, template, Ns) {
  values <- rowSums((component[, , 1] - component[, , 2]) ^ 2)
  values <- rep(values, each = Ns)
  vertices <- rbind(t(as.matrix(template[[2]])), 1)
  indices <- t(template[[3]])
  material <- values[1:dim(vertices)[2]]
  material <-
    as.integer(255 * (material - min(material)) / (max(material) - min(material)))
  paletteFunc <- colorRampPalette(c("blue", "green", "red"))
  
  palette <- paletteFunc(255)
  
  m <-
    tmesh3d(
      vertices = vertices,
      indices = indices,
      homogeneous = TRUE,
      material = list(palette[material])
    )
  shade3d(m, meshColor = "vertices")
}
library(tools)

plotAsymmetriesOnTemplate <- function(sym, template, Ns, saveDir="../results") {
  save <- getOption("rgl.useNULL")
  options(rgl.useNULL=TRUE)
  open3d();
  ids <- mfrow3d(1, 2)
  fa <- representComponent(sym$FA.component, template, Ns)
  title3d(main = "Fluctuating Asymmetry")
  next3d()
  da <- representComponent(sym$DA.component, template, Ns)
  title3d(main = "Directional Asymmetry")
  layout((1))
  widget <- rglwidget()
  if (interactive())
    widget
  
  # }
  # NOT RUN {
  # Save it to a file.  This requires pandoc
  
  dir.create(saveDir, showWarnings = FALSE, recursive = TRUE)
  filename <- file.path(saveDir, "interactive_asymmetries.html")
  htmlwidgets::saveWidget(rglwidget(), filename)
}


plotAsymmetriesOnTemplate(sym, template, Ns)
