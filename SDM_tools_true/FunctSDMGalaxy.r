#Rscript


##################################################################################################################################
####################### PAMPA Galaxy tools functions : Calculate metrics, compute GLM and plot   #################################
##################################################################################################################################


#### Modified by Coline ROYAUX for integrating within Galaxy-E

######################################### from virtualspecies package

synchroniseNA <- function(x) {
  if(canProcessInMemory(x, n = 2))
  {
    val <- getValues(x)
    if(any(is.na(val)))
    {
      NA.pos <- unique(which(is.na(val), arr.ind = T)[, 1])
    }
    val[NA.pos, ] <- NA
    x <- setValues(x, val)
    return(x)
  } else
  {
    x <- mask(x, calc(x, fun = sum))
    return(x)
  }
}

#########################################


