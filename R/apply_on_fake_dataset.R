require(geomorph)
library(stringr)

mydata <- read.csv("../SAMPLE_DATA/fake_dataset.csv", header = TRUE)

coords <- as.matrix(mydata[,2:dim(mydata)[2]])
sym_input.individual <- as.integer(sapply(strsplit(as.character(mydata$index),'_'), "[", 2))
sym_input.side <- as.integer(sapply(strsplit(as.character(mydata$index),'_'), "[", 3))
sym_input.rep <- as.integer(sapply(strsplit(as.character(mydata$index),'_'), "[", 4))

sym_input_df <- data.frame("individual"=sym_input.individual, "side"=sym_input.side, "rep"=sym_input.rep);
sym_coords <- aperm(array(coords, c(800, 3, 50)),c(3,2,1))

gpa <- gpagen(sym_coords, Proj=TRUE,PrinAxes = TRUE,max.iter=20)
summary(gpa)
sym <- bilat.symmetry(gpa,
                      ind = sym_input_df$individual,
                      side = sym_input_df$side,
                      replicate=sym_input_df$rep,
                      RRPP = FALSE,
                      iter = 50)
summary(sym)



head(toSave)

dim(mosquito$wingshape)




