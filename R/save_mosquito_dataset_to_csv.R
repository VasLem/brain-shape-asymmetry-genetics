require(geomorph)
data(mosquito)

toSave <- data.frame(array(aperm(mosquito$wingshape, c(3,2,1)), c(40, 36)))
names(toSave) <- rep(c("x","y"),18)
toSave$index <- paste("id", mosquito$ind, mosquito$side, mosquito$replicate, sep="_")
toSave <- toSave[,c(37,1:36)]
write.csv(toSave, file="../SAMPLE_DATA/mosquito.csv",row.names = FALSE,quote=FALSE)
data(mosquito)
gdf <- geomorph.data.frame(wingshape = mosquito$wingshape, ind=mosquito$ind, side=mosquito$side,
                           replicate=mosquito$replicate)
mosquito.gpa <- gpagen(gdf$wingshape, Proj=TRUE,PrinAxes = TRUE)
summary(mosquito.sym)

mosquito.sym <- bilat.symmetry(mosquito.gpa, ind = mosquito$ind, side = mosquito$side,
                               replicate = mosquito$replicate, object.sym = FALSE, RRPP = TRUE, iter = 50)
summary(mosquito.sym)
