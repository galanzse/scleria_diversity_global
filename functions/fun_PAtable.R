
# create PA.user.table for PA.strategy=user.defined in biomod2::BIOMOD_FormatingData

PAtable <- function(npres=NULL, PA.nb.rep=NULL, PA.nb.absences=NULL) {
  
  # define matrix
  PAtable <- matrix(nrow=npres+PA.nb.absences*PA.nb.rep, ncol=PA.nb.rep) %>% as.data.frame()

  # presence data, same for every run
  PAtable[1:npres,1:3] <- TRUE
  
  # pseudoabsences
  PAtable[c(1+npres):c(PA.nb.absences+npres),1] <- TRUE
  PAtable[c(PA.nb.absences+npres+1):c(PA.nb.absences*2+npres),2] <- TRUE
  PAtable[c(PA.nb.absences*2+npres+1):c(PA.nb.absences*3+npres),3] <- TRUE
  
  # fill NAs
  PAtable[is.na(PAtable)] <- FALSE
  
  return(PAtable)

}
