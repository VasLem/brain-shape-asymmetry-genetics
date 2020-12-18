require(geomorph)
data("mosquito")
library(Morpho)
# 18 2 40
data(mosquito)

wingshape <- mosquito$wingshape

library(tensor)
mosquito.gpa <- gpagen( mosquito$wingshape, Proj=TRUE,PrinAxes = TRUE)

avgShape <- mosquito.gpa$consensus

library(abind)


wingshape <- abind(mosquito$wingshape, array(0,c(18,1, 40)),along=2)                          
                                             
rotate <- function(s,theta=pi/3){
  p1 <- array(wingshape[s,,3])
  p2 <- array(wingshape[s,,1])
  p0 <- array(wingshape[s,,14])
  axisP1 <- p1 -  ( (p1 + p2) /2 - p0)
  
  axisP2 <-p2 -  ( (p1 + p2) /2 - p0)
  
  
  rotated <- t(rotaxis3d(t(wingshape[s,,]), pt1=axisP1, pt2=axisP2, theta=theta))
  rotated
}


rotated = array(0, c(18,3,40))
for (i in 1:18){
  rotated[i,,] = rotate(i)
}



#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

plot <- function(s)
{  
  rgl_init()
x <- rotated[s,1,1:18]
y <- rotated[s,2,1:18]
z <- rotated[s,3,1:18]
rgl.spheres(x,y,z, r = 0.002, color = "yellow")
rgl.bbox(color = "#333377") # Add bounding box decoration
x <- wingshape[s,1,1:18]
y <- wingshape[s,2,1:18]
z <- wingshape[s,3,1:18]
rgl.spheres(x,y,z, r = 0.002, color = "red")
}

# Make a scatter plot
plot(1)
ret = abind(wingshape, rotated,along=1)
toSave <- data.frame(array(aperm(ret, c(3,2,1)),c(40,18*3*2)))
names(toSave) <- rep(c("x","y","z"),18*2)
toSave$index <- paste("id", mosquito$ind, mosquito$side, mosquito$replicate, sep="_")
toSave <- toSave[,c(109,1:108)]
write.csv(toSave, file="../SAMPLE_DATA/mosquito3d.csv",row.names = FALSE,quote=FALSE)

# With the whole dataset
mdf <- data.frame(ind=mosquito$ind, side=mosquito$side,
                           replicate=mosquito$replicate)
mosquito.gpa <- gpagen(ret, Proj=TRUE,PrinAxes = TRUE)
summary(mosquito.gpa)

mosquito.sym <- bilat.symmetry(mosquito.gpa, ind = mdf$ind, side = mdf$side,
                               replicate = mdf$replicate, object.sym = FALSE, RRPP = TRUE, iter = 100)
summary(mosquito.sym)

# With the half dataset
mdf.half <- mdf[1:20,]
mosquito.gpaH <- gpagen(ret[,,1:20], Proj=TRUE,PrinAxes = TRUE)
summary(mosquito.gpaH)
mosquito.symH <- bilat.symmetry(mosquito.gpaH, ind = mdf.half$ind, side = mdf.half$side,
                               replicate = mdf.half$replicate, object.sym = FALSE, RRPP = TRUE, iter = 100)
summary(mosquito.symH)

# With double the dataset
mdf.double <- rbind(mdf,mdf)
mosquito.gpaD <- gpagen(abind(ret,ret,along=3), Proj=TRUE,PrinAxes = TRUE)
summary(mosquito.gpa)

mosquito.symD <- bilat.symmetry(mosquito.gpaD, ind = mdf.double$ind, side = mdf.double$side,
                               replicate = mdf.double$replicate, object.sym = FALSE, RRPP = TRUE, iter = 100)
summary(mosquito.sym)

