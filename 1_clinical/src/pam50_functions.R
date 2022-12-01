# Functions for PAM50 project

# loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
}

# define "not in" operator 
'%!in%' <- function(x,y)!('%in%'(x,y))

# 