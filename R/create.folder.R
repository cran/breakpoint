create.folder <-
function(i){
  x <- paste(getwd(),"/", i, sep="")
  dir.create(x)
  setwd(x)
}
