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

######################################### from biomod2 package v3.4.13

##' @title Convert species' probability of occurrence into binary presence-absence
##'   data using a predefined threshold
##'
##' @description Function that converts an object containing probability values
##'   into a binary presence-absence object according to a pre-defined threshold(s).
##'
##' @param data numeric vector, a \code{matrix}, a \code{data.frame}, a
##'   \code{RasterLayer} or a \code{RasterStack} containing the data to be
##'   converted
##' @param threshold numeric value or a vector containing the threshold to be
##'   used for converting data.
##'
##' @details
##'   If data is a vector or a raster object, then the threshold should be a
##'   numeric value. If data is matrix,dataframe or rasterStack, then the threshold
##'   should have, in theory, as many values as the number of columns or layers
##'   to transform.
##'   In the particular case that the data to convert is a \code{matrix}/\code{data.frame}
##'   with several columns or a \code{RasterStack} with several layers and the
##'   threshold is a single numeric value, the same threshold will be applied
##'   to all columns (resp. layers).
##'
##'
##' @return An object of the same class than \code{data} with binary (0 or 1) values,
##'   usually presence-absence.
##'
##' @author Wilfried Thuiller, Damien Georges
##'
##' @examples
##'   xx <- rnorm(50,10)
##'   yy <- BinaryTransformation(xx, 10)
##'
##'   cbind(xx,yy)
##'
##' @keywords models
##'
##' @export
##' @docType methods
##' @rdname BinaryTransformation-methods
setGeneric("BinaryTransformation",
           function(data, threshold){
             standardGeneric("BinaryTransformation")
           })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, data.frame-method
setMethod('BinaryTransformation', signature(data='data.frame'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='data.frame')")
  	FUN2 <- function(x,y){
  		moa <- apply((x>y),2,as.integer)
  		if(ncol(moa)==1) return(moa[,1])
  		else return(moa)
  	}
    if(is.numeric(threshold)){
      return(sweep(data.matrix(data), 2, threshold, FUN2))
    } else { ## return NAs
      return( matrix(NA, ncol=ncol(data), nrow=nrow(data), dimnames=dimnames(data)) )
    }

  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, matrix-method
setMethod('BinaryTransformation', signature(data='matrix'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='matrix')")
    data <- as.data.frame(data)
    return(BinaryTransformation(data, threshold))
  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, numeric-method
setMethod('BinaryTransformation', signature(data='numeric'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='numeric')")
    data <- as.data.frame(data)
    return(BinaryTransformation(data, threshold))
  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, array-method
setMethod('BinaryTransformation', signature(data='array'),
          function(data, threshold)
          {
            # cat("\n*** in setMethod('BinaryTransformation', signature(data='array')")
            if(length(dim(data)) == length(dim(threshold))){
              if(sum( dim(data)[-1] != dim(threshold)[-1] ) > 0 ){
                stop("data and threshold dimensions mismatch")
              }
            } else{
              if(sum( dim(data)[-1] != dim(threshold) ) > 0 ){
                stop("data and threshold dimensions mismatch")
              }
            }

            return(sweep(data,2:length(dim(data)),threshold,
                         function(x,y) {
                           if(!any(is.na(x))){
                             return(x>y)
                            } else {
                             return(rep(NA,length(x)) )}
                           } ))
          })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, RasterLayer-method
setMethod('BinaryTransformation', signature(data='RasterLayer'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='RasterLayer')")
    if(!is.na(threshold)){
      return(reclassify(data,c(-Inf,threshold,0, threshold,+Inf,1)))
    } else{ ## return a empty map (NA everywhere)
      return(reclassify(data,c(-Inf,Inf,NA)))
    }

  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, RasterStack-method
setMethod('BinaryTransformation', signature(data='RasterStack'),
  function(data, threshold)
  {
    if(length(threshold) == 1){
      threshold <- rep(threshold, raster::nlayers(data))
    }
    return(calc(data, function(x){x >= threshold}))
#     ## old version
#     StkTmp <- raster::stack()
#     for(i in 1:raster::nlayers(data)){
#       StkTmp <- raster::addLayer(StkTmp, BinaryTransformation(raster::subset(data,i,drop=TRUE), threshold[i]))
#     }
#     names(StkTmp) <- names(data)
#     return(StkTmp)
  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, RasterBrick-method
setMethod('BinaryTransformation', signature(data='RasterBrick'),
  function(data, threshold)
  {
    data <- raster::stack(data)
    return(BinaryTransformation(data, threshold))
  })
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Definition of models build with biomod2 to make it easy to access
# plot, .predict...
# Damien G. - 20/11/2012
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# NOTE -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This script will induce strong changes in biomod2. The version
# of biomod2 using objects defined there will be biomod2_2.x.x
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Formal Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('biomod2_model',
         representation(model_name = 'character',
                        model_class = 'character',
                        model_options = 'list',
                        model = 'ANY',
                        scaling_model = 'ANY',
                        resp_name = 'character',
                        expl_var_names = 'character',
                        expl_var_type = 'character',
                        expl_var_range = 'list',
                        model_evaluation = 'matrix',
                        model_variables_importance = 'matrix'),
         prototype = list(model_name = 'mySpecies_DataSet_RunName_myModelClass',
                   model_class = 'myModelClass',
                   model_options = list(),
                   model = list(),
                   scaling_model = list(),
                   resp_name = 'mySpecies',
                   expl_var_names = 'myRespVar',
                   expl_var_type = 'unknown',
                   expl_var_range = list(),
                   model_evaluation = matrix(),
                   model_variables_importance = matrix()),
         validity = function(object){

           # check that scaler is a glm if it is defined
           if(length(object@scaling_model))
             if(!inherits(object@scaling_model, c("glm", "lm"))) 
               return(FALSE)

           return(TRUE)
           } )

setMethod('show', signature('biomod2_model'),
          function(object){
            .bmCat("'biomod2_model'")
            cat("\n\t model name :", object@model_name, fill=.Options$width)
            cat("\n\t model class :", object@model_class, fill=.Options$width)
            cat("\n\t This model", ifelse(length(object@scaling_model), "has", "doesn't have"),"its own scaler", fill=.Options$width)

            cat("\n")
            cat("\n\t response modelled :", object@resp_name, fill=.Options$width)
            cat("\n\n\t explanatory variables used:", fill=.Options$width)
            cat("\n\t", "name", "\t", "type", "\t", "range", fill=.Options$width)
            for(i in 1: length(object@expl_var_names)){
              cat("\n\t", object@expl_var_names[i],"\t", object@expl_var_type[i], "\t", object@expl_var_range[[i]], fill=.Options$width)
            }

            cat("\n")
            cat("\n\t NOTE : ")
            cat("\n\t\t You can access 'formal' model with get_formal_model function")
            cat(ifelse(length(object@scaling_model), "\n\t\t You can access scaling model with get_scaling_model function\n", "\n"))

            .bmCat()
          })

# Models getters -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
setGeneric( "get_formal_model",
            def = function(object){
              standardGeneric( "get_formal_model" )
            } )

setMethod('get_formal_model', signature('biomod2_model'),
          function(object){
            return(object@model)
          })

setGeneric( "get_scaling_model",
            def = function(object){
              standardGeneric( "get_scaling_model" )
            } )

setMethod('get_scaling_model', signature('biomod2_model'),
          function(object){
            return(object@scaling_model)
          })

# Fuction to get variables ranges
get_var_type <- function(data){
    return(sapply(data,class))
}

get_var_range <- function(data){
  get_range <- function(x){
    if(is.numeric(x)){
      return(c(min=min(x,na.rm=T), max=max(x,na.rm=T)))
    }
    if(is.factor(x)){
      return(levels(x))
    }
  }
  xx <- lapply(data,get_range)
  names(xx) <- names(data)
  return(xx)
}

# Function to check new data range compatibility with calibrating data #
check_data_range <- function(model, new_data){
  ## TODO : remettre en marche cette fonction

#   # get calibration data caracteristics
#   expl_var_names <- model@expl_var_names
#   expl_var_type <- model@expl_var_type
#   expl_var_range <- model@expl_var_range
#
#   if(inherits(new_data, "Raster")){ ## raster data case =-=-=-=-=-=-=- #
#     # check var names compatibility
#     nd_expl_var_names <- names(new_data)
#     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
#       stop("calibration and projections variables names mismatch")
#     }
#     # reorder the stack
#     new_data <- raster::subset(new_data,expl_var_names)
#     # check var types compatibility (factors)
#     expl_var_fact <- (expl_var_type=='factor')
#     nd_expl_var_fact <- is.factor(new_data)
#     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
#       stop("calibration and projections variables class mismatch")
#     }
#     # check var range compatibility
#     ### remove all new factors
#     if(sum(expl_var_fact)>0){ ## there are factorial variables
#       for(fact_var_id in which(expl_var_fact)){
#         ## check if new factors occurs
#         nd_levels <- levels(raster::subset(new_data,fact_var_id))[[1]]
#         nd_levels <- as.character(nd_levels[,ncol(nd_levels)])
#         names(nd_levels) <- levels(raster::subset(new_data,fact_var_id))[[1]]$ID
#         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
#
#         ## detect new levels
#         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
#
#         if(length(new_levels)){
#           for(n_l in new_levels){
#             # remove points where out of range factors have been detected
#             new_data[subset(new_data,fact_var_id)[]==as.numeric(names(nd_levels)[which(nd_levels==n_l)])] <- NA
#           }
#           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
#         }
#       }
#     }
#     ## convert data to be sure to get RasterStack output
# #     new_data <- stack(new_data)
#   } else{ ## table data case -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#     # check var names compatibility
#     nd_expl_var_names <- colnames(new_data)
#     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
#       stop("calibration and projections variables names mismatch")
#     }
#     # reorder the stack
#     new_data <- new_data[,expl_var_names, drop=F]
#     # check var types compatibility (factors)
#     expl_var_fact <- (expl_var_type=='factor')
#     nd_expl_var_fact <- sapply(new_data,is.factor)
#
#     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
#       stop("calibration and projections variables class mismatch")
#     }
#     # check var range compatibility
#     ### remove all new factors
#     if(sum(expl_var_fact)>0){ ## there are factorial variables
#       for(fact_var_id in which(expl_var_fact)){
#         ## check if new factors occurs
#         nd_levels <- levels(new_data[,fact_var_id])
#         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
#
#         ## detect new levels
#         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
#
#         if(length(new_levels)){
#           # remove points where out of range factors have been detected
# #           new_data <- new_data[- which(new_data[,fact_var_id] %in% new_levels),]
#           new_data[which(new_data[,fact_var_id] %in% new_levels),] <- NA
#           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
#         }
#       }
#     }
#   }

  return(new_data)
}


# ANN Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setClass('ANN_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'ANN'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "nnet")) return(FALSE)
           return(TRUE)
           })

setMethod('predict', signature(object = 'ANN_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.ANN_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
               return(.predict.ANN_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.ANN_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

#   set.seed(555)
  proj <- predict(newdata, get_formal_model(object), type="raw")

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }

  return(proj)
}

.predict.ANN_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE



  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }


  set.seed(555)
  proj <- as.numeric( predict(get_formal_model(object), newdata[not_na_rows,,drop=F], type="raw") )

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# CTA Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('CTA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'CTA'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "rpart")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'CTA_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.CTA_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.CTA_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.CTA_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  set.seed(123) # to be able to refind our trees MAY BE BAD
  proj <- predict(newdata, model=get_formal_model(object), type='prob', index=2)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.CTA_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  set.seed(123)
  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]),type="prob")[,2])

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# FDA Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('FDA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'FDA'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "fda")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'FDA_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.FDA_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.FDA_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.FDA_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- predict(newdata, model=get_formal_model(object), type='posterior', index=2)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.FDA_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000

  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]),type = "posterior")[,2])

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# GAM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GAM_biomod2_model',
         representation(model_subclass='character'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GAM',
                   model_subclass = 'GAM_mgcv'),
         validity = function(object){
           # check model class
           if(! (object@model_subclass %in% c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv') ) )
             return(FALSE)

           if(object@model_subclass %in% c('GAM_mgcv','GAM_gam'))
             if(!inherits(object@model, c("gam")))
               return(FALSE)

           if(object@model_subclass == 'BAM_mgcv')
             if(!inherits(object@model, c("bam")))
               return(FALSE)

           return(TRUE)
         })

setMethod('predict', signature(object = 'GAM_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE


            if(object@model_subclass %in% c("GAM_mgcv","BAM_mgcv")){
              # cat("\n*** unloading gam package / loading mgcv package")
              if(isNamespaceLoaded("gam")){unloadNamespace("gam")}
              if(!isNamespaceLoaded("mgcv")){requireNamespace("mgcv", quietly = TRUE)}
            }

            if(object@model_subclass == "GAM_gam"){
              # cat("\n*** unloading mgcv package / loading gam package")
              if(isNamespaceLoaded("mgcv")){
                if(isNamespaceLoaded("caret")){unloadNamespace("caret")} ## need to unload caret before car
                if(isNamespaceLoaded("car")){unloadNamespace("car")} ## need to unload car before mgcv
                unloadNamespace("mgcv")
              }
              if(!isNamespaceLoaded("gam")){requireNamespace("gam", quietly = TRUE)}

            }

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.GAM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.GAM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.GAM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.GAM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows,,drop=FALSE])))
#   proj <- as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows,,drop=FALSE])))

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# GBM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GBM_biomod2_model',
         representation(n.trees_optim = 'numeric'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GBM',
                   n.trees_optim = 1000),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "gbm")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'GBM_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.GBM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.GBM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.GBM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

#   proj <- predict(newdata, model=get_formal_model(object), n.trees = object@n.trees_optim, type = "response")
  proj <- predict(newdata, model=get_formal_model(object), fun=gbm::predict.gbm, n.trees = object@n.trees_optim, type = "response")

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.GBM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), n.trees = object@n.trees_optim, type = "response"))

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# GLM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GLM_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GLM'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "glm")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'GLM_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.GLM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              ## transform matrix into dataframe to be able to use predict.glm fct
              if (inherits(newdata, 'matrix')) newdata <- as.data.frame(newdata)
              return(.predict.GLM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.GLM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.GLM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata[not_na_rows,,drop=FALSE]))

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# MARS Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('MARS_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MARS'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, c('earth', 'MARS', 'mars'))) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'MARS_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.MARS_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.MARS_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.MARS_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  ##' @note we have to handle separatly rasterstack depending on the presence or not
  ##' of factorial variable.
  fact.var <- is.factor(newdata)
  if(any(fact.var)){
    ## get factor levels
    fact.var.levels <- subset(levels(newdata), fact.var)
    proj <- calc(newdata, function(x) {
      xx <- data.frame(x)
      ## ensure that the data.frame has the right set of levels
      for(i in which(fact.var)){
        xx[[i]] <- factor(xx[[i]], levels = unlist(fact.var.levels[[i]]))
      }
      ## do the projection
      proj.out <- as.numeric(predict(get_formal_model(object), xx, type = 'response'))
      return(proj.out)
    })
  } else {
    proj <- predict(newdata, model = get_formal_model(object), type = 'response')
  }



  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.MARS_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), type = 'response'))


  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# MAXENT.Phillips Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('MAXENT.Phillips_biomod2_model',
         representation(model_output_dir = 'character'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MAXENT.Phillips'),
         validity = function(object){
           # check model class
#            if(sum(! ( c("randomForest.formula", "randomForest") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'MAXENT.Phillips_biomod2_model'),
          function(object, newdata, silent=TRUE, ...){
            args <- list(...)
            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.MAXENT.Phillips_biomod2_model.RasterStack(object, newdata, silent=TRUE, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.MAXENT.Phillips_biomod2_model.data.frame(object, newdata, silent=TRUE, ... ))
            } else{ stop("invalid newdata input") }

          })


.predict.MAXENT.Phillips_biomod2_model.RasterStack <- function(object, newdata, silent=TRUE,  ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  rm_tmp_files <- args$rm_tmp_files
  temp_workdir <- args$temp_workdir
  split.proj <- args$split.proj

#   if (is.null(temp_workdir)) temp_workdir <- paste("maxentWDtmp", format(Sys.time(), "%s"), sep="")
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if(is.null(split.proj)) split.proj <- 1

  MWD <- .Prepare.Maxent.Proj.WorkDir(Data = newdata, species.name = object@resp_name,
                                      silent = TRUE, split.proj = split.proj )

  # checking maxent.jar is present
  path_to_maxent.jar <- file.path(object@model_options$path_to_maxent.jar, "maxent.jar")
  if(!file.exists(path_to_maxent.jar)){
    path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
  }

  if(!silent) cat("\n\t\tRunning Maxent...")

  for(spl in 1:split.proj){
    maxent.command <- paste0("java ", ifelse(is.null(object@model_options$memory_allocated),"",paste0("-mx", object@model_options$memory_allocated, "m")),
                            " -cp ", "\"", path_to_maxent.jar, "\"",
                            " density.Project ",
                            "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                            "\"", MWD$m_workdir[[spl]], "\" ",
                            "\"", file.path(MWD$m_workdir[[spl]], "projMaxent.asc"), "\" ",
                            " doclamp=false visible=false autorun nowarnings notooltips")
    system(command = maxent.command, wait = TRUE, intern=TRUE)
  }

  if(!silent) cat("\n\t\tReading Maxent outputs...")

  ## get the list of projections part by part
  # check crs is not NA
  if(!is.na(projection(newdata))){
    proj.list <- lapply(file.path(unlist(MWD$m_workdir),"projMaxent.asc"), raster, RAT = FALSE, crs=projection(newdata))
  } else {
    proj.list <- lapply(file.path(unlist(MWD$m_workdir),"projMaxent.asc"), raster, RAT = FALSE)
  }
  ## merge all parts in a single raster
  if(length(proj.list) > 1){
    proj <- do.call(raster::merge, proj.list)
  } else {
    proj <- proj.list[[1]]
  }

#   # keep the coordinates ref system of new data
#   # TO DO => do it in .Prepare.Maxent.Proj.WorkDir()
#   proj.ref <- projection(newdata)

  #

  #   if(length(get_scaling_model(object))){
  #     names(proj) <- "pred"
  #     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  #   }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  } else if(!inMemory(proj)){
    proj <- readAll(proj) # to prevent from tmp files removing
  }




  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
#       unlink(x=file.path(object@resp_name, temp_workdir), recursive=TRUE, force=TRUE )
      .Delete.Maxent.WorkDir(MWD, silent=silent)
    }
  }

  return(proj)
}

.predict.MAXENT.Phillips_biomod2_model.data.frame <- function(object, newdata, silent=TRUE, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  temp_workdir <- args$temp_workdir
  rm_tmp_files <- args$rm_tmp_files
  xy <- args$xy

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
#   if (is.null(temp_workdir)) temp_workdir <- paste("maxentWDtmp", format(Sys.time(), "%s"), sep="")
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE

#   if( is.null(xy) ){
#     if( sum(c('x','y') %in% colnames(newdata) ) == 2 ){
#       coor_col <- c( which(colnames(newdata) == 'x'), which(colnames(newdata) == 'y') )
#       xy <- newdata[,coor_col]
#       newdata <- newdata[,- coor_col]
#     } else {
#       xy <- data.frame(x=rep(0,nrow(newdata)), y=rep(0,nrow(newdata)))
#     }
#   }

  ## no xy needed for models projections
  xy <- NULL

  ## check if na occurs in newdata cause they are not well supported
  not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})

  MWD <- .Prepare.Maxent.Proj.WorkDir(Data = as.data.frame(newdata[not_na_rows,,drop=FALSE]), xy = xy , species.name = object@resp_name, silent=T)
#   .Prepare.Maxent.Proj.WorkDir(Data = newdata, species.name = object@resp_name, proj.name = temp_workdir )

  # checking maxent.jar is present
  path_to_maxent.jar <- file.path(object@model_options$path_to_maxent.jar, "maxent.jar")
  if(!file.exists(path_to_maxent.jar)){
    path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
  }

  if(!silent) cat("\n\t\tRunning Maxent...")

  maxent.command <- paste0("java ", ifelse(is.null(object@model_options$memory_allocated),"",paste0("-mx", object@model_options$memory_allocated, "m")),
                           " -cp ", "\"", path_to_maxent.jar, "\"",
                           " density.Project ",
                           "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                           "\"", file.path(MWD$m_workdir, "Pred_swd.csv"), "\" ",
                           "\"", file.path(MWD$m_workdir, "projMaxent.asc") , "\" ",
                           "doclamp=false visible=false autorun nowarnings notooltips")

  system(command = maxent.command, wait = TRUE, intern=TRUE)


  if(!silent) cat("\n\t\tReading Maxent outputs...")
  proj <- as.numeric(read.asciigrid(file.path(MWD$m_workdir, "projMaxent.asc"))@data[,1])

  ## add original NAs in table
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      .Delete.Maxent.WorkDir(MWD, silent=silent)
    }
  }

#   if(length(get_scaling_model(object))){
#     proj <- data.frame(pred = proj)
#     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)

}

## MAXENT old Class ------------------------------------------------------------

##' @note Here is define a conversion of the MAXENT models until the biomod 2.xx
##'   This model class was equivqlent to the current 'MAXENT.Phillips_biomod2_model'

setClass('MAXENT_biomod2_model',
         contains = 'MAXENT.Phillips_biomod2_model',
         prototype = list(model_class = 'MAXENT'),
         validity = function(object){
           return(TRUE)
         })

## End MAXENT old Class --------------------------------------------------------

## Maxent Phillips 2 end ===========================================================================
setClass(
  'MAXENT.Phillips.2_biomod2_model',
  representation(),
  contains = 'biomod2_model',
  prototype = list(model_class = 'MAXENT.Phillips.2'),
  validity = function(object) checkmate::test_class(object@model, 'maxnet')
)

setMethod(
  'predict', 
  signature(object = 'MAXENT.Phillips.2_biomod2_model'),
  function(object, newdata, ...){
    args <- list(...)
    do_check <- args$do_check
    if(is.null(do_check)) do_check <- TRUE
    
    ## data checking
    if(do_check){
      newdata <- check_data_range(model=object, new_data=newdata)
    }
    
    if(inherits(newdata, 'Raster')){
      return(.predict.MAXENT.Phillips.2_biomod2_model.RasterStack(object, newdata, ... ))
    } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
      return(.predict.MAXENT.Phillips.2_biomod2_model.data.frame(object, newdata, ... ))
    } else { stop("invalid newdata input") }
    
  })

.predict.MAXENT.Phillips.2_biomod2_model.RasterStack <- function(object, newdata, ...){
  newdata.df <- 
    newdata %>%
    raster::as.matrix()
  
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- 
    predict(
      object = get_formal_model(object), 
      newdata = newdata.df,
      clamp = FALSE,
      type = 'logistic'
    )[, 1]
  
  if(length(get_scaling_model(object))){
    proj.to.scale <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj.to.scale)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  ## convert back to raster file
  proj.ras <- raster(newdata)
  proj.ras[apply(newdata.df, 1, function(.x) all(!is.na(.x)))] <- proj
  proj <- proj.ras
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }
  
  return(proj)
}

.predict.MAXENT.Phillips.2_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE
  
  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }
  
  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), type='logistic')[, 1])
  
  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }
  
  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}


## Maxent Phillips 2 end ===========================================================================

## MAXENT.Tsuruoka Class -------------------------------------------------------

# setClass('MAXENT.Tsuruoka_biomod2_model',
#          representation(),
#          contains = 'biomod2_model',
#          prototype = list(model_class = 'MAXENT.Tsuruoka'),
#          validity = function(object){
#            # check model class
#            if(sum(! ( c("maxent") %in% class(object@model) ) ) > 0) return(FALSE)
#            return(TRUE)
#          })

# setMethod('predict', signature(object = 'MAXENT.Tsuruoka_biomod2_model'),
#           function(object, newdata, ...){
#             args <- list(...)
# 
#             do_check <- args$do_check
#             if(is.null(do_check)) do_check <- TRUE
# 
#             ## data checking
#             if(do_check){
#               newdata <- check_data_range(model=object, new_data=newdata)
#             }
# 
#             if(inherits(newdata, 'Raster')){
#               return(.predict.MAXENT.Tsuruoka_biomod2_model.RasterStack(object, newdata, ... ))
#             } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
#               return(.predict.MAXENT.Tsuruoka_biomod2_model.data.frame(object, newdata, ... ))
#             } else{ stop("invalid newdata input") }
# 
#           })

# .predict.MAXENT.Tsuruoka_biomod2_model.RasterStack <- function(object, newdata, ...){
#   args <- list(...)
#   filename <- args$filename
#   overwrite <- args$overwrite
#   on_0_1000 <- args$on_0_1000
# 
#   if (is.null(overwrite)) overwrite <- TRUE
#   if (is.null(on_0_1000)) on_0_1000 <- FALSE
# 
#   proj <- calc(newdata, function(x) {
#     proj.out <- rep(NA, nrow(x))
#     x.no.na <- na.omit(x)
#     if(nrow(x.no.na)){
#       proj.not.na <- as.numeric(predict.maxent(get_formal_model(object), x.no.na)[, '1'])
#       proj.out[-attr(x.no.na, "na.action")] <- proj.not.na
#     }
#     return(proj.out)
#     })
# 
#   if(length(get_scaling_model(object))){
#     names(proj) <- "pred"
#     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }
# 
#   if(on_0_1000) proj <- round(proj*1000)
# 
#   # save raster on hard drive ?
#   if(!is.null(filename)){
#     cat("\n\t\tWriting projection on hard drive...")
#     if(on_0_1000){ ## projections are stored as positive integer
#       writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
#     } else { ## keep default data format for saved raster
#       writeRaster(proj, filename=filename, overwrite=overwrite)
#     }
#     proj <- raster(filename,RAT=FALSE)
#   }
# 
#   return(proj)
# }

# .predict.MAXENT.Tsuruoka_biomod2_model.data.frame <- function(object, newdata, ...){
#   args <- list(...)
#   on_0_1000 <- args$on_0_1000
#   omit.na <- args$omit.na
# 
#   if (is.null(on_0_1000)) on_0_1000 <- FALSE
#   if (is.null(omit.na)) omit.na <- FALSE
# 
#   ## check if na occurs in newdata cause they are not well supported
#   if(omit.na){
#     not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
#   } else {
#     not_na_rows <- rep(T, nrow(newdata))
#   }
# 
#   proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]))[,'1'])
# 
#   ## add original NAs in table if it needed
#   if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
#     tmp <- rep(NA,length(not_na_rows))
#     tmp[not_na_rows] <- proj
#     proj <- tmp
#     rm('tmp')
#   }
# 
#   if(length(get_scaling_model(object))){
#     proj <- data.frame(pred = proj)
#     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }
# 
#   if(on_0_1000) proj <- round(proj*1000)
# 
#   return(proj)
# }

## End MAXENT.Tsuruoka Class ---------------------------------------------------

# RF Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('RF_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'RF'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "randomForest")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'RF_biomod2_model'),
          function(object, newdata, ...){
            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.RF_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.RF_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.RF_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- predict(newdata, model=get_formal_model(object), type='prob', index=2)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }

  return(proj)
}

.predict.RF_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), type='prob')[,'1'])

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# SRE Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('SRE_biomod2_model',
         representation(extremal_conditions='data.frame'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'SRE'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'SRE_biomod2_model'),
          function(object, newdata, ...){
            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.SRE_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.SRE_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.SRE_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .sre.projection(NewData=newdata, ExtremCond=object@extremal_conditions)

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.SRE_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000

  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .sre.projection(NewData=newdata, ExtremCond=object@extremal_conditions)

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}

























# EM parent Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('biomod2_ensemble_model',
         representation(modeling.id='character'), ##maybe some additional args should be added here
         contains = 'biomod2_model',
         prototype = list(model_class = 'EM'),
         validity = function(object){
           return(TRUE)
         })

# EMmean Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMmean_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmean'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMmean_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){
            args <- list(...)

            do_check <- args$do_check

            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#            ## check if models are formal loaded
#             if(is.character(object@model)){
#               cat("\n*** models that will be loaded", object@model)
#               mem.change <- mem_change(model_tmp <- lapply(object@model, function(x){
#                 cat("\n***\t\t loading ", file.path(object@resp_name, "models", object@modeling.id, x))
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               }))
#               cat("\n*** mem change:", mem.change)
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMmean_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMmean_biomod2_model.data.frame(object, newdata, formal_predictions,  ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMmean_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
          function(mod, resp_name, modeling.id){
            ## check if model is loaded on memory
            if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
            return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
          }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  out <- calc(formal_predictions,
             function(x){
               m <- mean(x)
               if(on_0_1000) m <- round(m)
               return(m)
             },
             filename = filename, overwrite = TRUE)

  return(out)

}

.predict.EMmean_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                   function(mod, resp_name, modeling.id){
                                     ## check if model is loaded on memory
                                     if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                     return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                   } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  out <- rowMeans(formal_predictions, na.rm=T)

  if (on_0_1000){
    out <- round(out)
  }

  return(out)

}







# EMmedian Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMmedian_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMmedian_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMmedian_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMmedian_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMmedian_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  out <- calc(formal_predictions,
              function(x){
                m <- median(x)
                if(on_0_1000) m <- round(m)
                return(m)
              },
              filename = filename, overwrite = TRUE)

  return(out)

}

.predict.EMmedian_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod, resp_name, modeling.id){
                                   ## check if model is loaded on memory
                                   if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                   return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                 } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  out <- apply(formal_predictions, 1, median, na.rm=T)

  if (on_0_1000){
    out <- round(out)
  }

  return(out)

}








# EMcv Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMcv_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMcv_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)
            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMcv_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMcv_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMcv_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions
  mean_prediction <- args$mean_prediction

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  out <- calc(formal_predictions, cv,
              filename = filename, overwrite = TRUE, na.rm=TRUE, aszero=TRUE)

  return(out)

}

.predict.EMcv_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions
#   mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod, resp_name, modeling.id){
                                   ## check if model is loaded on memory
                                   if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                   return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                 } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

#   if(is.null(mean_prediction)){
#     # calculate mean of predictions
#     mean_prediction <- round(rowMeans(formal_predictions, na.rm=T))
#   }
#   # transforming 0 into Inf to produce null cv where mean is null
#   mean_prediction[mean_prediction==0] <- Inf
#
#   # calculate cv of formal models predictions
#   sd_prediction <- apply(formal_predictions,1,sd, na.rm=T)
#
#   return(round(sd_prediction / mean_prediction, 2))

  out <- apply(formal_predictions, 1, cv, na.rm=T, aszero=T)

#   if (on_0_1000){
#     out <- round(out)
#   }

  return(out)

}








# EMci Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMci_biomod2_model',
         representation(alpha = 'numeric',
                        side = 'character'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMci',
                   alpha = 0.05,
                   side = 'superior'),
         validity = function(object){
           if(!(object@side %in% c('inferior','superior'))) stop("side arg should be 'inferior' or 'superior")
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMci_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMci_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMci_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMci_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions
  mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving
  sd_prediction <- args$sd_prediction # mean of predictions should be given for time saving

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  if(is.null(mean_prediction)){
    mean_prediction <- calc(formal_predictions, mean)
  }

  if(is.null(sd_prediction)){
    sd_prediction <- calc(formal_predictions, sd)
  }

  ci_prediction <- switch(object@side,
                          inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                          superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ))

  # reclassify prediction to prevent from out of bounds prediction
  if (on_0_1000){
    ci_prediction <- reclassify(round(ci_prediction), c(-Inf,0,0, 1000,Inf,1000))
  } else {
    ci_prediction <- reclassify(ci_prediction, c(-Inf,0,0, 1,Inf,1))
  }

  return(ci_prediction)
}

.predict.EMci_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions
  mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving
  sd_prediction <- args$sd_prediction # mean of predictions should be given for time saving

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod, resp_name, modeling.id){
                                   ## check if model is loaded on memory
                                   if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                   return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                 } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  if(is.null(mean_prediction)){
    # calculate mean of predictions
    mean_prediction <- round(rowMeans(formal_predictions, na.rm=T))
  }

  if(is.null(sd_prediction)){
    # calculate cv of formal models predictions
    sd_prediction <- apply(formal_predictions,1,sd, na.rm=T)
  }

  ci_prediction <- switch(object@side,
                          inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                          superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ))

  # reclassify prediction to prevent from out of bounds prediction
  if (on_0_1000){
  ci_prediction <- round(ci_prediction)
  ci_prediction[ci_prediction > 1000] <- 1000
  ci_prediction[ci_prediction < 0] <- 0
  } else {
    ci_prediction[ci_prediction > 1] <- 1
    ci_prediction[ci_prediction < 0] <- 0
  }

  return(ci_prediction)
}







# EMca Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMca_biomod2_model',
         representation(tresholds = 'numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMca'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMca_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMca_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMca_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMca_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  if (on_0_1000){
    thresh <- object@tresholds
  } else { thresh <- object@tresholds / 1000 }

#   out <- raster::mean(BinaryTransformation(formal_predictions, thresh), na.rm=T)
  out <- calc(BinaryTransformation(formal_predictions, thresh),
                function(x){
                  m <- mean(x)
                  if(on_0_1000) m <- round(m * 1000)
                  return(m)
                  },
              filename = filename, overwrite = TRUE)

  return(out)

}

.predict.EMca_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod, resp_name, modeling.id){
                                   ## check if model is loaded on memory
                                   if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                   return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                 } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  if (on_0_1000){
    thresh <- object@tresholds
  } else { thresh <- object@tresholds / 1000 }

  out <- apply(as.data.frame(BinaryTransformation(formal_predictions, thresh)), 1, mean, na.rm=T)

  if (on_0_1000){
    out <- round(out * 1000)
  }

  return(out)

}

# EMwmean Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMwmean_biomod2_model',
         representation(penalization_scores='numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMwmean'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMwmean_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMwmean_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMwmean_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMwmean_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  filename <- args$filename
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

#   out <- sum(formal_predictions * object@penalization_scores)
  out <- calc(formal_predictions,
              function(x){
                wm <- sum(x * object@penalization_scores)
                if(on_0_1000) wm <- round(wm)
                return(wm)
              },
              filename = filename, overwrite = TRUE)

  return(out)

}

.predict.EMwmean_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod, resp_name, modeling.id){
                                   ## check if model is loaded on memory
                                   if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                   return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                 } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  out <- as.vector(as.matrix(formal_predictions) %*% object@penalization_scores)

  if (on_0_1000){
    out <- round(out)
  }

  return(out)

}

.Biomod.Models.loop <- function(X,
                                modeling.id,
                                Model,
                                Options,
                                VarImport,
                                mod.eval.method,
                                SavePred,
                                xy=NULL,
                                scal.models = TRUE){
  cat("\n\n-=-=-=- Run : ",X$name, '\n')
  res.sp.run <- list()

  for(i in 1:ncol(X$calibLines)){ # loop on RunEval
    cat('\n\n-=-=-=--=-=-=-',paste(X$name,dimnames(X$calibLines)[[2]][i],sep=""),'\n')

    res.sp.run[[dimnames(X$calibLines)[[2]][i]]] <- lapply(Model, .Biomod.Models,
                                                      Data = X$dataBM,
                                                      Options = Options,
                                                      calibLines = na.omit(X$calibLines[,i,]), ## transform 3D calibLines obj into a 1D vector
                                                      Yweights = na.omit(X$Yweights),
                                                      nam = paste(X$name,dimnames(X$calibLines)[[2]][i], sep=""),
                                                      VarImport = VarImport,
                                                      mod.eval.method = mod.eval.method,
                                                      evalData = X$evalDataBM,
                                                      SavePred = T,#SavePred,
                                                      xy = X$xy,
                                                      eval.xy = X$eval.xy,
                                                      scal.models = scal.models,
                                                      modeling.id = modeling.id)

    names(res.sp.run[[dimnames(X$calibLines)[[2]][i]]]) <- Model

  }

  return(res.sp.run)
}


.Biomod.Models <- function (Model, Data, Options, calibLines, Yweights, nam, VarImport = 0,
                            mod.eval.method = c('ROC','TSS','KAPPA'), evalData = NULL,
                            SavePred = FALSE,
                            xy = NULL, eval.xy = NULL, scal.models = TRUE, modeling.id = ''){

  ################################################################################################
  # 1. Print model running and getting model options =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  # check and get modified args if nececary
  args <- .Biomod.Models.check(Model, Data, Options, calibLines, Yweights, mod.eval.method, evalData, scal.models)

  if(is.null(args)){ # trouble in input data -> Not Run
    return(0)
  } else {
    Data <- args$Data
    Yweights <- args$Yweights
    evalLines <- args$evalLines
    Type <- args$Type
    criteria <- args$criteria
    Prev <- args$Prev
    mod.eval.method <- args$mod.eval.method
    evalData <- args$evalData
    scal.models <- args$scal.models
    resp_name <- args$resp_name
    expl_var_names <- args$expl_var_names
    compress.arg <- TRUE # ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }

  categorial_var <- unlist(sapply(expl_var_names, function(x){if(is.factor(Data[,x])) return(x) else return(NULL)} ))

  model_name <- paste(nam,'_',Model,sep="")


  # defining the function outputs
  ListOut <- list(evaluation = NULL,
                  var.import = NULL,
                  pred = NULL,
                  pred.eval = NULL,
                  calib.failure = NULL)


  # CTA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "CTA") {

    # converting cost argument
    if(is.null(Options@CTA$cost)){
      cost.tmp <- rep(1,(ncol(Data)-2))
    } else{
      cost.tmp <- Options@CTA$cost
    }
    if(Options@CTA$parms == 'default'){
      model.sp <- try( rpart(makeFormula(colnames(Data)[1],
                                         head(Data[,-c(1,ncol(Data)), drop=FALSE]),
                                         'simple', 0),
                             data = Data[calibLines,],
                             weights = Yweights,
                             method = Options@CTA$method,
                             cost = cost.tmp,
                             control = eval(Options@CTA$control)) )
    } else{
      model.sp <- try( rpart(makeFormula(colnames(Data)[1],
                                         head(Data[,-c(1,ncol(Data)), drop=FALSE]),
                                         'simple', 0),
                             data = Data[calibLines,],
                             weights = Yweights,
                             method = Options@CTA$method,
                             parms = Options@CTA$parms,
                             cost = cost.tmp,
                             control = eval(Options@CTA$control)) )
    }



    if( !inherits(model.sp,"try-error") ){
      # select best trees --------------- May be done otherway
      tr <- as.data.frame(model.sp$cptable)
      tr$xsum <- tr$xerror + tr$xstd
      tr <- tr[tr$nsplit > 0, ]
      Cp <- tr[tr$xsum == min(tr$xsum), "CP"]

      model.sp <- prune(model.sp, cp = Cp[length(Cp)])

      # creation of biomod2 model object
      model.bm <- new("CTA_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'CTA',
                      model_options = Options@CTA,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end CTA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #




  # GAM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "GAM"){

    # NOTE : To be able to take into account GAM options and weights we have to do a eval(parse(...))
    # it's due to GAM implementation ( using of match.call() troubles)

#     cat("\n\tUser defined control args building..")
#     user.control.list <- Options@GAM$control
#
#     if(Options@GAM$algo == 'GAM_gam'){
#       default.control.list <- gam::gam.control()
#     } else{
#       default.control.list <- mgcv::gam.control()
#     }
#
#     control.list <- lapply(names(default.control.list), function(x){
#       if(x %in% names(user.control.list)){
#         return(user.control.list[[x]])
#       } else {
#         return(default.control.list[[x]])
#       }
#     })
#     names(control.list) <- names(default.control.list)

    ### Old version
    if(Options@GAM$algo == 'GAM_gam'){ ## gam package
      # package loading
      if(isNamespaceLoaded("mgcv")){
        if(isNamespaceLoaded("caret")){unloadNamespace("caret")} ## need to unload caret before car
        if(isNamespaceLoaded("car")){unloadNamespace("car")} ## need to unload car before mgcv
        unloadNamespace("mgcv")
      }
      # if(!isNamespaceLoaded("gam")){requireNamespace("gam", quietly = TRUE)}
      requireNamespace("gam", quietly = TRUE)

      cat('\n\t> GAM (gam) modelling...')

      gamStart <- eval(parse(text=paste("gam::gam(",colnames(Data)[1] ,"~1 ," ,
                                        " data = Data[calibLines,,drop=FALSE], family = ", Options@GAM$family$family,"(link = '",Options@GAM$family$link,"')",#eval(Options@GAM$family),
                                        ", weights = Yweights[calibLines])" ,sep="")))
      model.sp <- try( gam::step.Gam(gamStart, .scope(Data[1:3,-c(1,ncol(Data))], "gam::s", Options@GAM$k),
                                     data = Data[calibLines,,drop=FALSE],
      #                                keep = .functionkeep,
                                     direction = "both",
                                     trace= Options@GAM$control$trace,
                                     control = Options@GAM$control))#eval(control.list)) )
    } else { ## mgcv package
      # package loading
#       if( ("package:gam" %in% search()) ){ detach("package:gam", unload=TRUE)}
#       if( ! ("package:mgcv" %in% search()) ){ require("mgcv",quietly=TRUE) }
      if(isNamespaceLoaded("gam")){unloadNamespace("gam")}
      if(!isNamespaceLoaded("mgcv")){requireNamespace("mgcv", quietly = TRUE)}


      if(is.null(Options@GAM$myFormula)){
        cat("\n\tAutomatic formula generation...")
        gam.formula <- makeFormula(resp_name,head(Data[,expl_var_names,drop=FALSE]),Options@GAM$type, Options@GAM$interaction.level, k=Options@GAM$k)
      } else{
        gam.formula <- Options@GAM$myFormula
      }

      if (Options@GAM$algo == 'GAM_mgcv'){
        cat('\n\t> GAM (mgcv) modelling...')
        model.sp <- try(mgcv::gam(gam.formula,
                                   data= Data[calibLines,,drop=FALSE],
                                   family= Options@GAM$family,
                                   weights = Yweights,
                                   control = Options@GAM$control))

      } else if (Options@GAM$algo == 'BAM_mgcv'){ ## big data.frame gam version
        cat('\n\t> BAM (mgcv) modelling...')
        model.sp <- try(mgcv::bam(gam.formula,
                                   data=Data[calibLines,,drop=FALSE],
                                   family=Options@GAM$family,
                                   weights = Yweights,
                                   control = Options@GAM$control))
      }
    }


    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("GAM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GAM',
                      model_subclass = Options@GAM$algo,
                      model_options = Options@GAM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end GAM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # GBM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "GBM") {

    model.sp <- try(gbm(formula = makeFormula(colnames(Data)[1],head(Data)[,expl_var_names,drop=FALSE], 'simple',0),
                        data = Data[calibLines,,drop=FALSE],
                        distribution = Options@GBM$distribution,
                        var.monotone = rep(0, length = ncol(Data)-2), # -2 because of removing of sp and weights
                        weights = Yweights,
                        interaction.depth = Options@GBM$interaction.depth,
                        n.minobsinnode = Options@GBM$n.minobsinnode,
                        shrinkage = Options@GBM$shrinkage,
                        bag.fraction = Options@GBM$bag.fraction,
                        train.fraction = Options@GBM$train.fraction,
                        n.trees = Options@GBM$n.trees,
                        verbose = Options@GBM$verbose,
                        #class.stratify.cv = Options@GBM$class.stratify.cv,
                        cv.folds = Options@GBM$cv.folds, 
                        n.cores = Options@GBM$n.cores)) ## to prevent from parallel issues

    if( !inherits(model.sp,"try-error") ){
      best.iter <- try(gbm.perf(model.sp, method = Options@GBM$perf.method , plot.it = FALSE))

      model.bm <- new("GBM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GBM',
                      n.trees_optim = best.iter,
                      model_options = Options@GBM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))

    }
  }
  # end GBM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # GLM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "GLM"){

    ## build the most complete model formula
    if(is.null(Options@GLM$myFormula)){
      glm.formula <- makeFormula(colnames(Data)[1],head(Data),Options@GLM$type, Options@GLM$interaction.level)
    } else{
      glm.formula <- Options@GLM$myFormula
    }

    if(Options@GLM$test != 'none'){
      ## make the model selection
      glmStart <- glm(eval(parse(text=paste(colnames(Data)[1],"~1",sep=""))),
                      data = Data[calibLines,,drop=FALSE],
                      family = Options@GLM$family,
                      control = eval(Options@GLM$control),
                      weights = Yweights[calibLines],
                      mustart = rep(Options@GLM$mustart, sum(calibLines)),
                      model = TRUE)

      ## remove warnings
      warn <- options('warn')
      options(warn=-1)
      model.sp <- try( stepAIC(glmStart,
                               glm.formula,
                               data = Data[calibLines,,drop=FALSE],
                               direction = "both", trace = FALSE,
                               k = criteria,
                               weights = Yweights[calibLines],
                               steps = 10000,
                               mustart = rep(Options@GLM$mustart, sum(calibLines))) )

      ## reexec warnings
      options(warn)

    } else {
      ## keep the total model
      model.sp <- try( glm(glm.formula,
                           data = cbind(Data[calibLines,,drop=FALSE],matrix(Yweights[calibLines], ncol=1, dimnames=list(NULL, "Yweights"))),
                           family = Options@GLM$family,
                           control = eval(Options@GLM$control),
                           weights = Yweights,
#                            mustart = rep(Options@GLM$mustart, sum(calibLines)),
                           model = TRUE) )
    }

    if( !inherits(model.sp,"try-error") ){
      # print the selected formula
      cat("\n\tselected formula : ")
      print(model.sp$formula, useSource=FALSE)



      model.bm <- new("GLM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GLM',
                      model_options = Options@GLM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end GLM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



#   # MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   if (Model == "MARS"){
#     ## deal with nk argument
#     ## if not defined, it will be setted up to default mars value i.e max(21, 2 * ncol(x) + 1)
#     nk <- Options@MARS$nk
#     if(is.null(nk)){
#       nk <- max(21, 2 * length(expl_var_names) + 1)
#     }
#
#     model.sp <- try( mars(x = Data[calibLines,expl_var_names,drop=FALSE],
#                           y = Data[calibLines,1],
#                           degree = Options@MARS$degree,
#                           nk = nk,
#                           penalty = Options@MARS$penalty,
#                           thresh = Options@MARS$thresh,
#                           prune = Options@MARS$prune,
#                           w = Yweights[calibLines]) )
#
#     if( !inherits(model.sp,"try-error") ){
#
#       model.bm <- new("MARS_biomod2_model",
#                       model = model.sp,
#                       model_name = model_name,
#                       model_class = 'MARS',
#                       model_options = Options@MARS,
#                       resp_name = resp_name,
#                       expl_var_names = expl_var_names,
#                       expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
#                       expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
#     }
#   }
#   # end MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

  # MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "MARS"){

    ## build the most complete model formula
    if(is.null(Options@MARS$myFormula)){
      mars.formula <- makeFormula(colnames(Data)[1],head(Data)[, -ncol(Data), drop = FALSE],Options@MARS$type, Options@MARS$interaction.level)
    } else{
      mars.formula <- Options@MARS$myFormula
    }

    ## deal with nk argument
    ## if not defined, it will be setted up to default mars value i.e max(21, 2 * ncol(x) + 1)
    nk <- Options@MARS$nk
    if(is.null(nk)){
      # nk <- max(21, 2 * length(expl_var_names) + 1)
      nk <- min(200, max(20, 2 * length(expl_var_names))) + 1
    }

    model.sp <- try(earth(formula = mars.formula,
                          data = Data[calibLines, , drop=FALSE],
                          weights = Yweights,
                          glm = list(family = binomial),
                          ncross = 0,
                          keepxy = FALSE,
                          # degree = Options@MARS$degree,
                          pmethod = Options@MARS$pmethod,
                          nprune = Options@MARS$nprune,
                          nk = nk,
                          penalty = Options@MARS$penalty,
                          thresh = Options@MARS$thresh))

    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("MARS_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'MARS',
                      model_options = Options@MARS,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # FDA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "FDA") {
    model.sp <- try( do.call(fda,
                             c( list( formula = makeFormula(colnames(Data)[1],head(Data)[,expl_var_names,drop=FALSE], 'simple',0),
                                      data = Data[calibLines,,drop=FALSE],
                                      method = eval(parse(text=call(Options@FDA$method))),
                                      weights = Yweights[calibLines] ),
                                Options@FDA$add_args) ) )

    if( !inherits(model.sp,"try-error") ){

      model.bm <- new("FDA_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'FDA',
                      model_options = Options@FDA,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end FDA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # ANN models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "ANN") {
    size = Options@ANN$size
    decay = Options@ANN$decay

    if(is.null(size) | is.null(decay) | length(size)>1 | length(decay)>1 ){
      ## define the size and decay to test
      if(is.null(size)) size <- c(2, 4, 6, 8)
      if(is.null(decay)) decay <- c(0.001, 0.01, 0.05, 0.1)

      ## do cross validation test to find the optimal values of size and decay parameters (prevent from overfitting)
      CV_nnet <- 
        .CV.nnet(
          Input = Data[,expl_var_names,drop=FALSE],
          Target = Data[calibLines,1],
          size = size,
          decay = decay,
          maxit = Options@ANN$maxit,
          nbCV = Options@ANN$NbCV,
          W = Yweights[calibLines]
        )

      ## get the optimised parameters values
      decay <- CV_nnet[1, 2]
      size <- CV_nnet[1,1]
    }

    model.sp <- try(nnet(formula = makeFormula(resp_name,head(Data[,expl_var_names,drop=FALSE]), 'simple',0),
                         data = Data[calibLines,,drop=FALSE],
                         size = size,
                         rang = Options@ANN$rang,
                         decay = decay,
                         weights=Yweights,
                         maxit = Options@ANN$maxit,
                         trace = FALSE))

    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("ANN_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'ANN',
                      model_options = Options@ANN,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # ANN models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # RF models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "RF") {
    if(Options@RF$do.classif){
      # defining occurences as factor for doing classification and not regression in RF
      Data <- Data %>% mutate_at(resp_name, factor)
    }

    if(Options@RF$mtry == 'default'){
      model.sp <- try(randomForest(formula = makeFormula(resp_name,head(Data), 'simple',0),
                                   data = Data[calibLines,],
                                   ntree = Options@RF$ntree,
                                   #mtry = ifelse(Options@RF$ntree == 'default', round((ncol(Data)-1)/2), Options@RF$ntree ),
                                   importance = FALSE,
                                   norm.votes = TRUE,
                                   strata = factor(c(0,1)),
                                   nodesize = Options@RF$nodesize,
                                   maxnodes = Options@RF$maxnodes) )
    } else {
      model.sp <- try(randomForest(formula = makeFormula(resp_name,head(Data), 'simple',0),
                                   data = Data[calibLines,],
                                   ntree = Options@RF$ntree,
                                   mtry = Options@RF$mtry,
                                   importance = FALSE,
                                   norm.votes = TRUE,
                                   strata = factor(c(0,1)),
                                   nodesize = Options@RF$nodesize,
                                   maxnodes = Options@RF$maxnodes) )
    }


    if(Options@RF$do.classif){
      # canceling occurences class modifications
      Data <- Data %>% mutate_at(resp_name, function(.x) .x %>% as.character() %>% as.numeric())
    }

    if( !inherits(model.sp,"try-error") ){

      model.bm <- new("RF_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'RF',
                      model_options = Options@RF,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end RF models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #





  # SRE models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "SRE"){
    model.sp <- try(sre(Response = Data[calibLines,1],
                        Explanatory = Data[calibLines,expl_var_names,drop=FALSE],
                        NewData = NULL,
                        Quant = Options@SRE$quant,
                        return_extremcond=TRUE))

    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("SRE_biomod2_model",
                      extremal_conditions = model.sp,
                      model_name = model_name,
                      model_class = 'SRE',
                      model_options = Options@SRE,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end SRE models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #





  # MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "MAXENT.Phillips"){
    MWD <- .Prepare.Maxent.WorkDir(Data, xy, calibLines, nam, VarImport = 0,
									evalData, eval.xy, species.name=resp_name,
									modeling.id = modeling.id,
                                   	background_data_dir = Options@MAXENT.Phillips$background_data_dir)

    # run MaxEnt:
    cat("\n Running Maxent...")
    maxent.cmd <- paste0("java ",
                        ifelse(is.null(Options@MAXENT.Phillips$memory_allocated),"",paste("-mx",Options@MAXENT.Phillips$memory_allocated,"m",sep="")),
                        " -jar ", file.path(Options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"),
                        " environmentallayers=\"", MWD$m_backgroundFile,
                        "\" samplesfile=\"", MWD$m_speciesFile,
                        "\" projectionlayers=\"", gsub(", ",",",toString(MWD$m_predictFile)),
                        "\" outputdirectory=\"", MWD$m_outdir, "\"",
                        " outputformat=logistic ",
                        #                            "jackknife maximumiterations=",Options@MAXENT.Phillips$maximumiterations,
                        ifelse(length(categorial_var),
                               paste(" togglelayertype=",categorial_var, collapse=" ",sep=""),
                               ""),
                        " redoifexists",
                        " visible=", Options@MAXENT.Phillips$visible,
                        " linear=", Options@MAXENT.Phillips$linear,
                        " quadratic=", Options@MAXENT.Phillips$quadratic,
                        " product=", Options@MAXENT.Phillips$product,
                        " threshold=", Options@MAXENT.Phillips$threshold,
                        " hinge=", Options@MAXENT.Phillips$hinge,
                        " lq2lqptthreshold=", Options@MAXENT.Phillips$lq2lqptthreshold,
                        " l2lqthreshold=", Options@MAXENT.Phillips$l2lqthreshold,
                        " hingethreshold=", Options@MAXENT.Phillips$hingethreshold,
                        " beta_threshold=", Options@MAXENT.Phillips$beta_threshold,
                        " beta_categorical=", Options@MAXENT.Phillips$beta_categorical,
                        " beta_lqp=", Options@MAXENT.Phillips$beta_lqp,
                        " beta_hinge=", Options@MAXENT.Phillips$beta_hinge,
                        " betamultiplier=", Options@MAXENT.Phillips$betamultiplier,
                        " defaultprevalence=", Options@MAXENT.Phillips$defaultprevalence,
                        " autorun nowarnings notooltips noaddsamplestobackground")

    system(command = maxent.cmd, wait = TRUE, intern = TRUE,
           ignore.stdout = FALSE, ignore.stderr = FALSE)

    model.bm <- new("MAXENT.Phillips_biomod2_model",
                    model_output_dir = MWD$m_outdir,
                    model_name = model_name,
                    model_class = 'MAXENT.Phillips',
                    model_options = Options@MAXENT.Phillips,
                    resp_name = resp_name,
                    expl_var_names = expl_var_names,
                    expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                    expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))

    # for MAXENT.Phillips predicitons are calculated in the same time than models building to save time.
    cat("\n Getting predictions...")
    g.pred <- try(round(as.numeric(read.csv(MWD$m_outputFile)[,3]) * 1000))

    # remove tmp dir
    .Delete.Maxent.WorkDir(MWD)
  }
  # end MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  # MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(Model == "MAXENT.Phillips.2")
  {
    # browser()
    model.sp <- 
      try(
        maxnet::maxnet(
          p = Data %>% filter(calibLines) %>% pull(resp_name), 
          data = Data %>% filter(calibLines) %>% select_at(expl_var_names)
          # f = if(!is.null(Options@MAXENT.Phillips.2@))
        )
      )
    
    
    if( !inherits(model.sp,"try-error") )
    {
      model.bm <- 
        new(
          "MAXENT.Phillips.2_biomod2_model",
          model = model.sp,
          model_name = model_name,
          model_class = 'MAXENT.Phillips.2',
          model_options = Options@MAXENT.Phillips.2,
          resp_name = resp_name,
          expl_var_names = expl_var_names,
          expl_var_type = get_var_type(Data %>% filter(calibLines) %>% select_at(expl_var_names)),
          expl_var_range = get_var_range(Data %>% filter(calibLines) %>% select_at(expl_var_names))
        )
    }
  }
  # end MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-= #

  # # MAXENT.Tsuruoka models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  # if(Model == "MAXENT.Tsuruoka"){
  #   model.sp <- try(stop('MAXENT.Tsuruoka is depreacated(because maxent package is not maintained anymore)'))
  #   # model.sp <- try(maxent::maxent(feature_matrix = Data[calibLines, expl_var_names, drop = FALSE],
  #   #                                code_vector = as.factor(Data[calibLines, 1]),
  #   #                                l1_regularizer = Options@MAXENT.Tsuruoka$l1_regularizer,
  #   #                                l2_regularizer = Options@MAXENT.Tsuruoka$l2_regularizer,
  #   #                                use_sgd = Options@MAXENT.Tsuruoka$use_sgd,
  #   #                                set_heldout = Options@MAXENT.Tsuruoka$set_heldout,
  #   #                                verbose = Options@MAXENT.Tsuruoka$verbose))
  # 
  #   if( !inherits(model.sp,"try-error") ){
  #     model.bm <- new("MAXENT.Tsuruoka_biomod2_model",
  #                     model = model.sp,
  #                     model_name = model_name,
  #                     model_class = 'MAXENT.Tsuruoka',
  #                     model_options = Options@MAXENT.Tsuruoka,
  #                     resp_name = resp_name,
  #                     expl_var_names = expl_var_names,
  #                     expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
  #                     expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
  #   }
  # }
  # # end of MAXENT.Tsuruoka models creation -=-=-=-=-=-=-=-=-=-=-=- #

  # make prediction =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if((Model != "MAXENT.Phillips")){
    g.pred <- try(predict(model.bm, Data[, expl_var_names, drop = FALSE], on_0_1000 = TRUE))
  }


  # scale or not predictions =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(scal.models & !inherits(g.pred,'try-error')){
    cat("\n\tModel scaling...")
    #     model.bm@scaling_model <- try(.scaling_model(g.pred/1000, Data[, 1])) ## without weigths
    model.bm@scaling_model <- try( .scaling_model(g.pred/1000, Data[, 1, drop = TRUE], weights= Yweights )) ## with weights
    g.pred <- try(predict(model.bm, Data[,expl_var_names,drop=FALSE], on_0_1000=TRUE))
  }

  # check predictions existance and stop execution if not ok -=-=- #
  test_pred_ok <- TRUE
  if (inherits(g.pred,"try-error")) { # model calibration or prdiction failed
    test_pred_ok <- FALSE
    cat("\n*** inherits(g.pred,'try-error')")
  } else if (sum(!is.na(g.pred))<=1){ # only NA predicted
    test_pred_ok <- FALSE
    cat("\n*** only NA predicted")
  } else if(length(unique(na.omit(g.pred))) <=1){ # single value predicted
    test_pred_ok <- FALSE
    cat("\n*** single value predicted")
  }

  if(test_pred_ok){
    # keep the model name
    ListOut$ModelName <- model_name
  } else{
    # keep the name of uncompleted modelisations
    cat("\n   ! Note : ", model_name, "failed!\n")
    ListOut$calib.failure = model_name
    return(ListOut) ## end of function.
  }

  # make prediction on evaluation data =-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(!is.null(evalData)){
    g.pred.eval <- try(predict(model.bm, evalData[,expl_var_names,drop=FALSE], on_0_1000=TRUE))
  }

  # save predictions -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(SavePred){
    ListOut$pred <- g.pred
    if(exists("g.pred.eval"))
      ListOut$pred.eval <- g.pred.eval
  }


  # Model evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(length(mod.eval.method) > 0){
    cat("\n\tEvaluating Model stuff...")

    ## Check no NA in g.pred to avoid evaluation failures
    na_cell_id <- which(is.na(g.pred))
    if(length(na_cell_id)){
#       g.pred.without.na <- g.pred[-na_cell_id]
      evalLines <- evalLines[!(evalLines %in% na_cell_id)]
      cat('\n\tNote : some NA occurs in predictions')
    } #else {
#       g.pred.without.na <- g.pred
#     }

    cross.validation <- 
      sapply(
        mod.eval.method,
        function(.x){
          Find.Optim.Stat(
            Stat = .x,
            Fit = g.pred[evalLines],
            Obs = Data %>% filter(evalLines) %>% pull(1)
          )
        }
      )

    rownames(cross.validation) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")

    if(exists('g.pred.eval')){

      ## Check no NA in g.pred to avoid evaluation failures
      na_cell_id <- which(is.na(g.pred.eval))
      if(length(na_cell_id)){
        g.pred.eval.without.na <- g.pred.eval[-na_cell_id]
        evalData <- evalData[-na_cell_id,]
        cat('\n\tNote : some NA occurs in evaluation predictions')
      } else {
        g.pred.eval.without.na <- g.pred.eval
      }

      true.evaluation <- sapply(mod.eval.method,
                                function(x){
                                  return( Find.Optim.Stat(Stat = x,
                                                          Fit = g.pred.eval.without.na,
                                                          Obs = evalData[,1],
                                                          Fixed.thresh = cross.validation["Cutoff",x]) )
                                })


      cross.validation <- rbind(cross.validation["Testing.data",], true.evaluation)

      rownames(cross.validation) <- c("Testing.data","Evaluating.data","Cutoff","Sensitivity", "Specificity")
    }

    ListOut$evaluation <- t(round(cross.validation,digits=3))

    ## store results
    cross.validation <- t(round(cross.validation,digits=3))
    ListOut$evaluation <- cross.validation
    model.bm@model_evaluation <- cross.validation

    ## remove useless objects
    rm(list=c('cross.validation') )
  }
  # End evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #




  # Variables Importance -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (VarImport > 0){ # do Varimp stuff
    cat("\n\tEvaluating Predictor Contributions...", "\n")
    variables.importance <- variables_importance(model.bm, Data[, expl_var_names,drop=FALSE], nb_rand=VarImport)
    model.bm@model_variables_importance <- variables.importance$mat
    ## we stored only the mean of variables importance run
    ListOut$var.import <- round(rowMeans(variables.importance$mat, na.rm=T),digits=3)
    ## remove useless objects
    rm(list=c('variables.importance') )
  }
#   if (VarImport > 0){ # do Varimp stuff
#     # Create data.frame vhere corelation between predictions with and without permutation will be
#     # stored
#
#     cat("\n\tEvaluating Predictor Contributions...", "\n")
#     VarImpTable <- matrix(data = 0, nrow = VarImport, ncol = length(expl_var_names))
#     dimnames(VarImpTable) <- list(paste('rand', 1:VarImport, sep=""), expl_var_names)
#
#     for(vari in expl_var_names){
#       for (run in 1:VarImport) {
#         ## create a new dataset with interest variable suffled
#         TempDS <- Data[, expl_var_names,drop=FALSE]
#         TempDS[, vari] <- sample(TempDS[, vari])
#
#         if(Model != "MAXENT.Phillips"){
#           ## make projection on suffled dataset
#           shuffled.pred <- try(predict(model.bm, TempDS, on_0_1000=TRUE))
#         } else{
#           ## for MAXENT.Phillips, we have created all the permutation at model building step
#           shuffled.pred <- try(round(as.numeric(read.csv(file.path(model.bm@model_output_dir, paste(nam, vari, run, "swd.csv", sep="_")))[,3])*1000) )
#           ## scal suffled.pred if necessary
#           if(length(getScalingModel(model.bm))){
#             shuffled.pred <- try( round(.testnull(object = getScalingModel(model.bm), Prev = 0.5 , dat = data.frame(pred = shuffled.pred/1000) ) *1000) )
#             #               shuffled.pred <- round(as.numeric(predict(getScalingModel(model.bm), shuffled.pred/1000))*1000)
#           }
#           ## remove useless files on hard drive
#           file.remove(list.files(path=model.bm@model_output_dir,
#                                  pattern=paste(nam, vari, run, "swd", sep="_"),
#                                  full.names=TRUE))
#         }
#
#         ## test if differences exist between the 2 vectors
#         # check predictions existance and stop execution if not ok -=-=- #
#         test_shuffled.pred_ok <- TRUE
#         if (inherits(shuffled.pred,"try-error")) { # model calibration or prdiction failed
#           test_shuffled.pred_ok <- FALSE
#         } else if (sum(!is.na(shuffled.pred))<=1){ # only NA predicted
#           test_shuffled.pred_ok <- FALSE
#         } else if(length(unique(na.omit(shuffled.pred))) <=1){ # single value predicted
#           test_shuffled.pred_ok <- FALSE
#         } else if(length(shuffled.pred)!= length(g.pred)){
#           test_shuffled.pred_ok <- FALSE
#         }
#
#         if(!test_shuffled.pred_ok){
#           cat("\n   ! Note : ", model_name, "variable importance for",vari,run,"failed!\n")
#           VarImpTable[run,vari] <- 0
#         } else{
#           if(sum( g.pred != shuffled.pred, na.rm=T) == 0){
#             VarImpTable[run,vari] <- 0
#           } else {
#             ## calculate correlation between vectors as proxy for variables importance
#             VarImpTable[run,vari] <- 1 - max(round(cor(x=g.pred, y=shuffled.pred, use="pairwise.complete.obs", method="pearson"),digits=3),0,na.rm=T)
#           }
#         }
#
#       }
#     }
#
#     ## store results
#     model.bm@model_variables_importance <- VarImpTable
#     ## we stored only the mean of variables importance run
#     ListOut$var.import <- round(apply(VarImpTable, 2, mean, na.rm=T),digits=3)
#   }
  # End Variables Importance -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # Model saving step =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  assign(x=paste( nam, Model, sep = "_"),
         value= model.bm)
  save(list=paste( nam, Model, sep = "_"),
       file=file.path(resp_name, "models", modeling.id, paste( nam, Model, sep = "_")),
       compress=compress.arg)


  # End model saving step =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

  return(ListOut)
}


.Biomod.Models.check <- function(Model, Data, Options, calibLines, Yweights, mod.eval.method, evalData, scal.models, criteria=NULL, Prev=NULL){
  # get species and expanatory variables names
  resp_name <- colnames(Data)[1]
  expl_var_names <- colnames(Data)[-1]

  # replace Pseudo absences selected (NA) into true absences (0).. for model computing purpose

  if(sum(is.na(Data[,1])))
    Data[which(is.na(Data[,1])),1] <- 0

  # Calib Lines checking
  # & Test if there is absences AND presences in data given
  if (sum(!calibLines)>0){ # data are splited into 2 set one for evaluation and an other for evaluation stuff
    evalLines <- !calibLines
    if(sum(Data[calibLines,1] == 0 ) == 0 || sum(Data[calibLines,1] == 0 ) == sum(calibLines) ||
      sum(Data[evalLines,1] == 0) == 0 || sum(Data[evalLines,1] == 0) == sum(evalLines)){
      warning(paste(colnames(Data)[1], " ", Model," was switched off because of no both
                   presences and absences data given",sep=""), immediate.=T)
     return(NULL)
    }
  } else { # all values are taken for cali and valid stuff -----> Not so good but usefull for little data set
    evalLines <- calibLines
    if(sum(Data[,1] == 0 ) == 0 || sum(Data[,1] == 0 ) == nrow(Data)){
      warning(paste(colnames(Data)[1], " ", Model," was switched off because of no both
                   presences and absences data given (full model)",sep=""), immediate.=T)
      return(NULL)
    }
  }

  # weights checking
  if(is.null(Yweights)){
    Yweights <- rep(1,nrow(Data))
  }

  if(Model %in% c('GBM', 'CTA', 'ANN', 'FDA', 'GAM', 'MARS')){ # this models required data and weights to be in a same datdaset
    Data <- cbind(Data,Yweights)
  }

  # scaling parameter checking
  # never scal SRE
  if(Model == "SRE") scal.models <- FALSE
  # always scal ANN, FDA
  if(Model %in% c("ANN", "FDA") ) scal.models <- TRUE


  # models options checking and printing
  if (Model == "GLM"){
    cat('\nModel=GLM')
    if(!is.null(Options@GLM$myFormula)){
      cat('\n\tformula = ', paste(Options@GLM$myFormula[2],Options@GLM$myFormula[1],Options@GLM$myFormula[3]))
    } else{
      cat(' (',Options@GLM$type,'with', ifelse(Options@GLM$interaction.level == 0, 'no interaction )', paste('order',Options@GLM$interaction.level,'interaction level )')))
    }

    if(Options@GLM$test == "AIC"){
      criteria <- 2
      cat("\n\tStepwise procedure using AIC criteria")
    } else if(Options@GLM$test == "BIC"){
      criteria <- log(ncol(Data))
      cat("\n\tStepwise procedure using BIC criteria")
    } else if(Options@GLM$test == "none"){
      criteria <- 0
      cat("\n\tNo stepwise procedure")
      cat("\n\t! You might be confronted to models convergence issues !")
    }

  }

  if (Model == "GBM") {
    cat("\nModel=Generalised Boosting Regression \n")
    cat("\t", Options@GBM$n.trees, "maximum different trees and ", Options@GBM$cv.folds,
        " Fold Cross-Validation")
    set.seed(456) # to be able to refind our trees MAY BE BAD
  }

  if (Model == "GAM") {
    cat("\nModel=GAM")
    cat("\n\t",Options@GAM$algo,"algorithm chosen")
    #         cat("\t", Options@GAM$spline, " Degrees of smoothing")
    #         if(is.null(Yweights)) Yweights <- rep(1,nrow(Data))
    #         Data <- cbind(Data,Yweights)
  }

  if (Model == "CTA") {
    cat("\nModel=Classification tree \n")
    cat("\t", Options@CTA$control$xval, "Fold Cross-Validation")
    set.seed(123) # to be able to refind our trees MAY BE BAD
  }

  if (Model == "ANN") {
    cat("\nModel=Artificial Neural Network \n")
    cat("\t", Options@ANN$NbCV, "Fold Cross Validation + 3 Repetitions")
    #         cat("\tCalibration and evaluation phase: Nb of cross-validations: ",
    #             ncol(Ids), "\n")
    set.seed(555) # to be able to refind our trees MAY BE BAD
  }

  if (Model == "SRE")
    cat("\nModel=Surface Range Envelop")

  if (Model == "FDA"){
    cat("\nModel=Flexible Discriminant Analysis")
  }

  if (Model == "MARS"){
    cat("\nModel=Multiple Adaptive Regression Splines")
    if(!is.null(Options@MARS$myFormula)){
      cat('\n\tformula = ', paste(Options@MARS$myFormula[2],Options@MARS$myFormula[1],Options@MARS$myFormula[3]))
    } else{
      cat(' (',Options@MARS$type,'with', ifelse(Options@MARS$interaction.level == 0, 'no interaction )', paste('order',Options@MARS$interaction.level,'interaction level )')))
    }
    cat("\n")
  }

  if (Model == "RF"){
    cat("\nModel=Breiman and Cutler's random forests for classification and regression")
    set.seed(71)
  }

  if(Model == 'MAXENT.Phillips'){
    cat('\nModel=MAXENT.Phillips')
  }
  
  if(Model == 'MAXENT.Phillips.2'){
    cat('\nModel=MAXENT.Phillips (maxnet)')
  }

  # if(Model == 'MAXENT.Tsuruoka'){
  #   cat('\nModel=MAXENT.Tsuruoka')
  # }

  #     if (Model == "GLM" | Model == "GAM")
  #         Prev <- sum(DataBIOMOD[, i + Biomod.material$NbVar])/nrow(DataBIOMOD)
  ## not exactly same as before
  if (Model == "GLM" | Model == "GAM"){
    Prev <- sum(Data[,1], na.rm=T)/length(Data[,1])
  }

  # Evaluation Check
  available.eval.meth <- c('ROC','KAPPA','TSS','ACCURACY','BIAS','POD','FAR','POFD','SR','CSI',
                           'ETS','HK','HSS','OR','ORSS')

  #   if( Model %in% c('SRE') ) available.eval.meth <- available.eval.meth[which(available.eval.meth!='ROC')]
  if(sum(!(mod.eval.method %in% available.eval.meth)) > 0 ){
    warnings(paste(toString(mod.eval.method[!which(mod.eval.method %in% available.eval.meth)]),
                   ' were switched off !', sep='' ), imediate = TRUE)
  }
  mod.eval.method <- mod.eval.method[which(mod.eval.method %in% available.eval.meth)]

  return(list(Data=Data,
              Yweights=Yweights,
              evalLines=evalLines,
              criteria=criteria,
              Prev=Prev,
              mod.eval.method=mod.eval.method,
              evalData=evalData,
              scal.models=scal.models,
              resp_name=resp_name,
              expl_var_names=expl_var_names))

}
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# BIOMOD objects definition
# Damien Georges
# 09/02/2012
# v2e("sp", quietly=TRUE)
requireNamespace("raster", quietly=TRUE)
requireNamespace(".0
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This file defines the BIOMOD objects and all their methods
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# We choose here to create monospecific objects to make all procedures and parallelising easier
requireNamespacrasterVis", quietly=TRUE)

# 0. Generic Functions definition -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
setGeneric("get_predictions",
           function(obj, ...){
             standardGeneric("get_predictions")
           })

setGeneric("get_projected_models",
           function(obj, ...){
             standardGeneric("get_projected_models")
           })

setGeneric("get_evaluations",
           function(obj, ...){
             standardGeneric("get_evaluations")
           })

setGeneric("get_calib_lines",
           function(obj, ...){
             standardGeneric("get_calib_lines")
           })

setGeneric("get_variables_importance",
           function(obj, ...){
             standardGeneric("get_variables_importance")
           })

setGeneric("get_options",
           function(obj, ...){
             standardGeneric("get_options")
           })

setGeneric("get_formal_data",
           function(obj, ...){
             standardGeneric("get_formal_data")
           })

setGeneric("get_built_models",
           function(obj, ...){
             standardGeneric("get_built_models")
           })

setGeneric("get_needed_models",
           function(obj, ...){
             standardGeneric("get_needed_models")
           })

setGeneric("get_kept_models",
           function(obj, ...){
             standardGeneric("get_kept_models")
           })

setGeneric("load_stored_object",
           function(obj, ...){
             standardGeneric("load_stored_object")
           })

# setGeneric("RemoveProperly",
#            function(obj, obj.name=deparse(substitute(obj)), ...){
#              standardGeneric("RemoveProperly")
#            })

setGeneric("free",
           function(obj, ...){
             standardGeneric("free")
           })

setGeneric( ".Models.prepare.data",
            def = function(data, ...){
              standardGeneric( ".Models.prepare.data" )
            } )


# 1. The BIOMOD.formated.data -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this object is the basic one

# 1.1 Class Definition
setClass("BIOMOD.formated.data",
         representation(sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
                        #                         data.counting = "matrix",
                        data.env.var = "data.frame",
                        data.mask = "RasterStack",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame"),
         validity = function(object){ return(TRUE) } )

# 1.2 Constructors
# if( !isGeneric( "BIOMOD.formated.data" ) ) {
setGeneric( "BIOMOD.formated.data",
            def = function(sp, env, ...){
              standardGeneric( "BIOMOD.formated.data" )
            } )
# }

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='data.frame' ),
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE, data.mask=NULL ){
            if(is.null(data.mask)) data.mask <- raster::stack()

            if(is.null(eval.sp)){
              BFD <- new('BIOMOD.formated.data',
                         coord=xy,
                         data.species=sp,
                         data.env.var=env,
                         sp.name=sp.name,
                         data.mask=data.mask,
                         has.data.eval=FALSE)
            } else{
              BFDeval <- BIOMOD.formated.data(sp=eval.sp,
                                              env=eval.env,
                                              xy=eval.xy,
                                              sp.name=sp.name)

              if(raster::nlayers(BFDeval@data.mask)>0){
                data.mask.tmp <- try(raster::addLayer(data.mask,BFDeval@data.mask))
                if( !inherits(data.mask.tmp,"try-error")){
                  data.mask <- data.mask.tmp
                  names(data.mask) <- c("calibration","validation")
                }
              }

              BFD <- new('BIOMOD.formated.data',
                         coord=xy,
                         data.species=sp,
                         data.env.var=env,
                         sp.name=sp.name,
                         data.mask=data.mask,
                         has.data.eval=TRUE,
                         eval.coord = BFDeval@coord,
                         eval.data.species = BFDeval@data.species,
                         eval.data.env.var = BFDeval@data.env.var )


              rm('BFDeval')
            }
            if(na.rm){
              rowToRm <- unique(unlist(lapply(BFD@data.env.var,function(x){return(which(is.na(x)))})))
              if(length(rowToRm)){
                cat("\n\t\t\t! Some NAs have been automatically removed from your data")
                BFD@coord <- BFD@coord[-rowToRm,,drop=FALSE]
                BFD@data.species <- BFD@data.species[-rowToRm]
                BFD@data.env.var <- BFD@data.env.var[-rowToRm,,drop=FALSE]
              }
              if(BFD@has.data.eval){
                rowToRm <- unique(unlist(lapply(BFD@eval.data.env.var,function(x){return(which(is.na(x)))})))
                if(length(rowToRm)){
                  cat("\n\t\t\t! Some NAs have been automatically removed from your evaluation data")
                  BFD@eval.coord <- BFD@eval.coord[-rowToRm,,drop=FALSE]
                  BFD@eval.data.species <- BFD@eval.data.species[-rowToRm]
                  BFD@eval.data.env.var <- BFD@eval.data.env.var[-rowToRm,,drop=FALSE]
                }
              }
            }



            # count data occutances
            #     BFD@data.counting <- matrix(c(sum(BFD@data.species, na.rm=TRUE),sum(BFD@data.species==0, na.rm=TRUE)),
            #                             ncol=1,nrow=2, dimnames=list(c("nb.pres","nb.abs"),c("data.species") ) )
            #
            #     if(BFD@has.data.eval){
            #       BFD@data.counting <- cbind(BFD@data.counting,c(sum(BFD@eval.data.species, na.rm=TRUE),sum(BFD@eval.data.species==0, na.rm=TRUE)))
            #       colnames(BFD@data.counting)[ncol(BFD@data.counting)] <- "eval.data.species"
            #     }

            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp='data.frame'),
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
            if(ncol(sp) > 1 ){
              stop("Invalid response variable")
            }
            sp <- as.numeric(unlist(sp))
            BFD <- BIOMOD.formated.data(sp,env,xy,sp.name, eval.sp, eval.env, eval.xy, na.rm=na.rm)
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='matrix' ),
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
            env <- as.data.frame(env)
            BFD <- BIOMOD.formated.data(sp,env,xy,sp.name, eval.sp, eval.env, eval.xy, na.rm=na.rm)
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='RasterStack' ),
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
            categorial_var <- names(env)[raster::is.factor(env)]

            # take the same eval environemental variables than calibrating ones
            if(!is.null(eval.sp)){
              if(is.null(eval.env)){
                #         eval.env_levels <- levels(eval.env)
                eval.env <- as.data.frame(extract(env,eval.xy))
                if(length(categorial_var)){
                  for(cat_var in categorial_var){
                    eval.env[,cat_var] <- as.factor(eval.env[,cat_var])
                  }
                }
              }
            }

            if(is.null(xy)) xy <- as.data.frame(coordinates(env))

            data.mask = reclassify(raster::subset(env,1,drop=T), c(-Inf,Inf,-1))
            data.mask[cellFromXY(data.mask,xy[which(sp==1),])] <- 1
            data.mask[cellFromXY(data.mask,xy[which(sp==0),])] <- 0
            data.mask <- raster::stack(data.mask)
            names(data.mask) <- sp.name

            #     env_levels <- levels(env)
            env <- as.data.frame(extract(env,xy, factors=T))

            if(length(categorial_var)){
              for(cat_var in categorial_var){
                env[,cat_var] <- as.factor(env[,cat_var])
              }
            }

            BFD <- BIOMOD.formated.data(sp,env,xy,sp.name,eval.sp, eval.env, eval.xy, na.rm=na.rm, data.mask=data.mask)
            return(BFD)
          }
)

# 1.3 Other Functions
# if( !isGeneric( "plot" ) ) {
#   setGeneric( "plot",
#               def = function(x, ...){
#   	                  standardGeneric( "plot" )
#                       } )
# }

setMethod('plot', signature(x='BIOMOD.formated.data', y="missing"),
          function(x,coord=NULL,col=NULL){
            if(raster::nlayers(x@data.mask)>0){
              requireNamespace("rasterVis")

              ## check that there is some undefined areas to prevent from strange plotting issues
              if(min(cellStats(x@data.mask,min)) == -1){ # there is undifined area
                ## define the breaks of the color key
                my.at <- seq(-1.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(-1,1,by=1)
                ## define the labels
                my.lab <- c("undifined","absences","presences")
                ## define colors
                my.col.regions = c("lightgrey","red4","green4")
                ## defined cuts
                my.cuts <- 2
              } else{ # no undefined area.. remove it from plot
                ## define the breaks of the color key
                my.at <- seq(-0.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(0,1,by=1)
                ## define the labels
                my.lab <- c("absences","presences")
                ## define colors
                my.col.regions = c("red4","green4")
                ## defined cuts
                my.cuts <- 1
              }


              levelplot(x@data.mask, at=my.at, cuts=my.cuts, margin=T, col.regions=my.col.regions,
                        main=paste(x@sp.name,"datasets"),
                        colorkey=list(labels=list(
                          labels=my.lab,
                          at=my.labs.at)))

            } else{
              # coordinates checking
              if(is.null(coord)){
                if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }

              # colors checking
              if(is.null(col) | length(col) < 3){
                col = c('green', 'red', 'grey')
              }

              # plot data
              # all points (~mask)

              plot(x=x@coord[,1], y=x@coord[,2], col=col[3], xlab = 'X', ylab = 'Y',
                   main = paste(x@sp.name, sep=""), pch=20 )
              # presences
              points(x=x@coord[which(x@data.species == 1),1],
                     y=x@coord[which(x@data.species == 1),2],
                     col=col[1],pch=18)
              # true absences
              points(x=x@coord[which(x@data.species == 0),1],
                     y=x@coord[which(x@data.species == 0),2],
                     col=col[2],pch=18)

            }



          })

##' @rdname BIOMOD.formated.data-objects
##' @docType method
##' @aliases show, BIOMOD.formated.data-method
setMethod('show', signature('BIOMOD.formated.data'),
          function(object){
            .bmCat("'BIOMOD.formated.data'")
            cat("\nsp.name = ", object@sp.name, fill=.Options$width)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ',
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n', fill=.Options$width)
            print(summary(object@data.env.var))

            if(object@has.data.eval){
              cat("\n\nEvaluation data :", fill=.Options$width)
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                  sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ',
                  sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
              cat("\n\n", fill=.Options$width)
              print(summary(object@eval.data.env.var))
            }

            .bmCat()
          })

# 2. The BIOMOD.formated.data.PA =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this class inherits from BIOMOD.formated.data and have one more slot 'PA' giving PA selected

# 2.1 Class Definition
setClass("BIOMOD.formated.data.PA",
         contains = "BIOMOD.formated.data",
         representation(PA.strategy='character', PA = 'data.frame'),
         validity = function(object){
           return(TRUE)
         })

# 2.2 Constructors
# if( !isGeneric( "BIOMOD.formated.data.PA" ) ){
#   setGeneric( "BIOMOD.formated.data.PA",
#               def = function(sp, env, PA.NbRep, ...){
#                 standardGeneric( "BIOMOD.formated.data.PA" )
#               })
# }

# setMethod('BIOMOD.formated.data.PA',signature(sp='ANY',
#                                               env='ANY',
#                                               PA.NbRep='integer'),

BIOMOD.formated.data.PA <-  function(sp, env, xy, sp.name,
                                     eval.sp=NULL, eval.env=NULL, eval.xy=NULL,
                                     PA.NbRep=1,
                                     PA.strategy='random',
                                     PA.nb.absences = NULL,
                                     PA.dist.min = 0,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     PA.table = NULL,
                                     na.rm=TRUE){

  if(inherits(env,'Raster')){
    categorial_var <- names(env)[raster::is.factor(env)]
  }  else categorial_var <- NULL


  # take the same eval environemental variables than calibrating ones
  if(!is.null(eval.sp)){
    if(is.null(eval.env)){
      if(inherits(env,'Raster')){
        eval.env <- as.data.frame(extract(env,eval.xy))
        if(length(categorial_var)){
          for(cat_var in categorial_var){
            eval.env[,cat_var] <- as.factor(eval.env[,cat_var])
          }
        }
      } else{
        stop("No evaluation explanatory variable given")
      }
    }
  }

  # convert sp in spatial obj
  if(is.numeric(sp)){
    if(is.null(xy)){
      sp <- SpatialPointsDataFrame(matrix(0,ncol=2,nrow=length(sp)), data.frame(sp),match.ID=FALSE)
    } else{
      sp <- SpatialPointsDataFrame(data.matrix(xy), data.frame(sp),match.ID=FALSE)
    }

  }

  pa.data.tmp <- .pseudo.absences.sampling(sp = sp,
                                           env = env,
                                           nb.repet = PA.NbRep,
                                           strategy = PA.strategy,
                                           nb.points = PA.nb.absences,
                                           distMin = PA.dist.min,
                                           distMax = PA.dist.max,
                                           quant.SRE = PA.sre.quant,
                                           PA.table = PA.table)

  if(!is.null(pa.data.tmp)){

    if(length(categorial_var)){
      for(cat_var in categorial_var){
        pa.data.tmp$env[,cat_var] <- as.factor(pa.data.tmp$env[,cat_var])
      }
    }


    if(na.rm){
      rowToRm <- unique(unlist(lapply(pa.data.tmp$env,function(x){return(which(is.na(x)))})))
      if(length(rowToRm)){
        cat("\n\t\t\t! Some NAs have been automatically removed from your data")
        pa.data.tmp$xy <- pa.data.tmp$xy[-rowToRm,,drop=FALSE]
        pa.data.tmp$sp <- pa.data.tmp$sp[-rowToRm, drop=FALSE]
        pa.data.tmp$env <- pa.data.tmp$env[-rowToRm,,drop=FALSE]
        pa.data.tmp$pa.tab <- pa.data.tmp$pa.tab[-rowToRm,,drop=FALSE]
        #         cat("\n***\n")
        #         cat("\ndim(xy) <-", dim(pa.data.tmp$xy))
        #         cat("\nclass(sp) <- ", class(pa.data.tmp$sp))
        #         cat("\ndim(sp) <-", dim(pa.data.tmp$sp))
        #         cat("\ndim(env) <-", dim(pa.data.tmp$env))
        #         cat("\ndim(pa.tab) <-", dim(pa.data.tmp$pa.tab))

      }
    }


    # data counting
    #     pa.data.tmp$data.counting <- apply(pa.data.tmp$pa.tab,2,function(x){nbPres <- sum(pa.data.tmp$sp[x],na.rm=T) ; return(c(nbPres,sum(x)-nbPres))})
    #     colnames(pa.data.tmp$data.counting) <- colnames(pa.data.tmp$pa.tab)

    BFD <- BIOMOD.formated.data(sp=pa.data.tmp$sp,
                                env=pa.data.tmp$env,
                                xy=as.data.frame(pa.data.tmp$xy),
                                sp.name=sp.name,
                                eval.sp=eval.sp,
                                eval.env=eval.env,
                                eval.xy=eval.xy,
                                na.rm=na.rm) # because check is already done

    if(inherits(env,'Raster')){

      ## create data.mask for ploting
      data.mask.tmp <- reclassify(raster::subset(env,1), c(-Inf,Inf,-1))
      data.mask <- stack(data.mask.tmp)
      xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp==1), , drop = FALSE]
      xy_abs <- pa.data.tmp$xy[which(pa.data.tmp$sp==0), , drop = FALSE]
      if(nrow(xy_pres)){
        data.mask[cellFromXY(data.mask.tmp, xy_pres)] <- 1
      }
      if(nrow(xy_abs)){
        data.mask[cellFromXY(data.mask.tmp, xy_abs)] <- 0
      }
      names(data.mask) <- "input_data"

      ## add eval data
      if(BFD@has.data.eval){
        ### TO DO

      }

      for(pa in 1:ncol(as.data.frame(pa.data.tmp$pa.tab))){
        data.mask.tmp2 <- data.mask.tmp

        xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp==1 & as.data.frame(pa.data.tmp$pa.tab)[,pa]), , drop = FALSE]
        xy_abs <- pa.data.tmp$xy[which( (pa.data.tmp$sp!=1 | is.na(pa.data.tmp$sp)) & as.data.frame(pa.data.tmp$pa.tab)[,pa]), , drop = FALSE]

        if(nrow(xy_pres)){
          id_pres <- cellFromXY(data.mask.tmp, xy_pres)
          data.mask.tmp2[id_pres] <- 1
        }

        if(nrow(xy_abs)){
          id_abs <- cellFromXY(data.mask.tmp, xy_abs)
          data.mask.tmp2[id_abs] <- 0
        }

        data.mask <- addLayer(data.mask, data.mask.tmp2)
      }

      names(data.mask) <- c("input_data", colnames(as.data.frame(pa.data.tmp$pa.tab)))

    } else{
      data.mask <- stack()
    }

    BFDP <- new('BIOMOD.formated.data.PA',
                sp.name = BFD@sp.name,
                coord = BFD@coord,
                #                 data.counting = cbind(BFD@data.counting,pa.data.tmp$data.counting) ,
                data.env.var = BFD@data.env.var,
                data.species = BFD@data.species,
                data.mask = data.mask,
                has.data.eval = BFD@has.data.eval,
                eval.coord = BFD@eval.coord,
                eval.data.species = BFD@eval.data.species,
                eval.data.env.var = BFD@eval.data.env.var,
                PA = as.data.frame(pa.data.tmp$pa.tab),
                PA.strategy = PA.strategy)

    rm(list='BFD')
  } else {
    cat("\n   ! PA selection not done", fill=.Options$width)

    BFDP <- BIOMOD.formated.data(sp=sp,
                                 env=env,
                                 xy=xy,
                                 sp.name=sp.name,
                                 eval.sp=eval.sp,
                                 eval.env=eval.env,
                                 eval.xy=eval.xy)

  }

  rm(list = "pa.data.tmp" )

  return(BFDP)

}


# 2.3 other functions
setMethod('plot', signature(x='BIOMOD.formated.data.PA', y="missing"),
          function(x,coord=NULL,col=NULL){

            if(raster::nlayers(x@data.mask)>0){
              requireNamespace("rasterVis")

              ## check that there is some undefined areas to prevent from strange plotting issues
              if(min(cellStats(x@data.mask,min)) == -1){ # there is undifined area
                ## define the breaks of the color key
                my.at <- seq(-1.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(-1,1,by=1)
                ## define the labels
                my.lab <- c("undifined","absences","presences")
                ## define colors
                my.col.regions = c("lightgrey","red4","green4")
                ## defined cuts
                my.cuts <- 2
              } else{ # no undefined area.. remove it from plot
                ## define the breaks of the color key
                my.at <- seq(-0.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(0,1,by=1)
                ## define the labels
                my.lab <- c("absences","presences")
                ## define colors
                my.col.regions = c("red4","green4")
                ## defined cuts
                my.cuts <- 1
              }


              levelplot(x@data.mask, at=my.at, cuts=my.cuts, margin=T, col.regions=my.col.regions,
                        main=paste(x@sp.name,"datasets"),
                        colorkey=list(labels=list(
                          labels=my.lab,
                          at=my.labs.at)))

            } else{
              # coordinates checking
              if(is.null(coord)){
                if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }

              # colors checking
              if(is.null(col) | length(col) < 3){
                col = c('green', 'red', 'orange', 'grey')
              }

              # plot data
              par(mfrow=c(.CleverCut(ncol(x@PA)+1)))
              # all points (~mask)
              plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
                   main = paste(x@sp.name," original data", sep=""), pch=20 )
              # presences
              points(x=x@coord[which(x@data.species == 1),1],
                     y=x@coord[which(x@data.species == 1),2],
                     col=col[1],pch=18)
              # true absences
              points(x=x@coord[which(x@data.species == 0),1],
                     y=x@coord[which(x@data.species == 0),2],
                     col=col[2],pch=18)

              # PA data
              for(i in 1:ncol(x@PA)){
                # all points (~mask)
                plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
                     main = paste(x@sp.name," Pseudo Absences ", i, sep=""), pch=20 )
                # presences
                points(x=x@coord[(x@data.species == 1) & x@PA[,i],1],
                       y=x@coord[(x@data.species == 1) & x@PA[,i],2],
                       col=col[1],pch=18)
                # true absences
                points(x=x@coord[(x@data.species == 0) & x@PA[,i],1],
                       y=x@coord[(x@data.species == 0) & x@PA[,i],2],
                       col=col[2],pch=18)
                # PA
                points(x=x@coord[is.na(x@data.species) & x@PA[,i],1],
                       y=x@coord[is.na(x@data.species) & x@PA[,i],2],
                       col=col[3],pch=18)
              }
            } })

##' @rdname BIOMOD.formated.data.PA-objects
##' @docType method
##' @aliases show, BIOMOD.formated.data.PA-method
setMethod('show', signature('BIOMOD.formated.data.PA'),
          function(object){
            .bmCat("'BIOMOD.formated.data.PA'")
            cat("\nsp.name = ", object@sp.name,fill=.Options$width)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ',
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n', fill=.Options$width)
            print(summary(object@data.env.var))

            if(object@has.data.eval){
              cat("\n\nEvaluation data :", fill=.Options$width)
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                  sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ',
                  sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
              cat("\n\n", fill=.Options$width)
              print(summary(object@eval.data.env.var))
            }

            cat("\n\n", ncol(object@PA), 'Pseudo Absences dataset available (', colnames(object@PA),") with ",
                sum(object@PA[,1], na.rm=T) - sum(object@data.species, na.rm=TRUE), 'absences in each (true abs + pseudo abs)', fill=.Options$width)
            .bmCat()
          })


##' @rdname BIOMOD.Model.Options-objects
##' @name BIOMOD.Model.Options-class
##' @docType class
##' @aliases   BIOMOD.Model.Options-class
##' 
##' @title BIOMOD_ModelingOptions outputs objects class
##' 
##' @description
##' BIOMOD.Model.Options objects are created, used and returned
##' by BIOMOD functions. These objects will contains for each
##' model support within \pkg{biomod2}, a set of options that
##' users can change. Please refer to 
##' \code{\link[biomod2]{BIOMOD_ModelingOptions}} for further
##'  details. 
##'   
##' - output of: \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' - input of:  \code{\link[biomod2]{BIOMOD_Modeling}}
##' 
##' @param object init list of options
##' 
##' @details  
##' Please refer to \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' for each model arguments supported.
##' 
##' @slot GLM "list", list of GLM supported options
##' @slot GBM "list", list of GBM supported options 
##' @slot GAM "list", list of GAM supported options
##' @slot CTA "list", list of CTA supported options
##' @slot ANN "list", list of ANN supported options
##' @slot SRE "list", list of SRE supported options
##' @slot FDA "list", list of FDA supported options
##' @slot MARS "list", list of MARS supported options
##' @slot RF "list", list of RF supported options
##' @slot MAXENT.Phillips "list", list of MAXENT.Phillips
##'   supported options
##' @slot MAXENT.Phillips.2 "list", list of maxnet
##'   supported options
##'   
##' @author Damien Georges
##' @seealso \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' @keywords models
##' @keywords options
##' 
##' @examples
##' showClass("BIOMOD.Model.Options")
##' 
setClass(
  "BIOMOD.Model.Options",
  representation(
    GLM = "list",
    GBM = "list",
    GAM = "list",
    CTA = "list",
    ANN = "list",
    SRE = "list",
    FDA = "list",
    MARS = "list",
    RF = "list",
    MAXENT.Phillips = "list",
    MAXENT.Phillips.2 = "list"
  ),
  prototype(
    GLM = 
      list( 
        type = 'quadratic',
        interaction.level = 0,
        myFormula = NULL,
        test = 'AIC',
        family = binomial(link = 'logit'),
        mustart = 0.5,
        control = glm.control(maxit = 50)
      ),
    GBM = 
      list(
        distribution = 'bernoulli',
        n.trees = 2500,
        interaction.depth = 7,
        n.minobsinnode = 5,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        train.fraction = 1,
        cv.folds = 3,
        keep.data = FALSE,
        verbose = FALSE,
        # class.stratify.cv = 'bernoulli',
        perf.method = 'cv',
        n.cores = 1
      ),
    GAM = 
      list( 
        algo = "GAM_mgcv",
        type = "s_smoother",
        k = NULL,
        interaction.level = 0,
        myFormula = NULL,
        family = binomial(link = 'logit'),
        control = list(epsilon = 1e-06, trace = FALSE, maxit = 100),
        method = "GCV.Cp",
        optimizer = c("outer", "newton"),
        select = FALSE,
        knots = NULL,
        paraPen = NULL
      ),
    CTA = 
      list(
        method = 'class',
        parms = 'default',
        # control = rpart.control(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25),
        control = 
          list(
            xval = 5, 
            minbucket = 5, 
            minsplit = 5,
            cp = 0.001, 
            maxdepth = 25
          ),
        cost = NULL 
      ),
    ANN = 
      list(
        NbCV = 5,
        size = NULL,
        decay = NULL,
        rang = 0.1,
        maxit = 200
      ),
    SRE = 
      list(
        quant = 0.025
      ),
    FDA = 
      list(
        method = 'mars',
        add_args = NULL
      ),
    MARS = 
      list(
        type = 'simple',
        interaction.level = 0,
        myFormula = NULL,
        # degree = 1,
        nk = NULL,
        penalty = 2,
        thresh = 0.001,
        nprune = NULL,
        pmethod = 'backward'
      ),
    RF = 
      list(
        do.classif = TRUE,
        ntree = 500,
        mtry = 'default',
        nodesize = 5,
        maxnodes= NULL
      ),
    MAXENT.Phillips = 
      list(
        path_to_maxent.jar = getwd(),
        memory_allocated = 512,
        background_data_dir = 'default',
        maximumbackground = 'default',
        maximumiterations = 200,
        visible = FALSE,
        linear = TRUE,
        quadratic = TRUE,
        product = TRUE,
        threshold = TRUE,
        hinge = TRUE,
        lq2lqptthreshold = 80,
        l2lqthreshold = 10,
        hingethreshold = 15,
        beta_threshold = -1.0,
        beta_categorical = -1.0,
        beta_lqp = -1.0,
        beta_hinge = -1.0,
        betamultiplier = 1,
        defaultprevalence = 0.5
      ),
    MAXENT.Phillips.2 = 
      list(
        myFormula = NULL,
        regmult = 1,
        regfun = maxnet::maxnet.default.regularization
      )
  ),
  validity = 
    function(object){
      test <- TRUE
      ## GLM ##
      if(!(object@GLM$type %in% c('simple','quadratic','polynomial','user.defined'))){ cat("\nGLM$type must be 'simple',  'quadratic', 'polynomial' or 'user.defined'"); test <- FALSE}
      if(!is.numeric(object@GLM$interaction.level)){ cat("\nGLM$interaction.level must be a integer"); test <- FALSE } else{
             if(object@GLM$interaction.level < 0 | object@GLM$interaction.level%%1!=0){ cat("\nGLM$interaction.level must be a positive integer"); test <- FALSE }
           }
           if(!is.null(object@GLM$myFormula)) if(!inherits(object@GLM$myFormula, "formula")){ cat("\nGLM$myFormula must be NULL or a formula object"); test <- FALSE }
           if(!(object@GLM$test %in% c('AIC','BIC','none'))){ cat("\nGLM$test must be 'AIC','BIC' or 'none'"); test <- FALSE}
           fam <- 'none'
           if(!inherits(object@GLM$family, "family")){ cat("\nGLM$family must be a valid family object"); test <- FALSE }
           if(!is.list(object@GLM$control)){cat("\nGLM$control must be a list like that returned by glm.control"); test <- FALSE}


           ## GBM ##
           if(! object@GBM$distribution %in% c("bernoulli","huberized", "multinomial", "adaboost")){cat("\nGBM$distribution must be 'bernoulli', 'huberized', 'multinomial'or  'adaboost'") }

           if(!is.numeric(object@GBM$n.trees)){ cat("\nGBM$n.trees must be a integer"); test <- FALSE } else{
             if(object@GBM$n.trees < 0 | floor(object@GBM$n.trees) != object@GBM$n.trees){ cat("\nGBM$n.trees must be a positive integer"); test <- FALSE }
           }

           if(!is.numeric(object@GBM$interaction.depth)){ cat("\nGBM$interaction.depth must be a integer"); test <- FALSE } else{
             if(object@GBM$interaction.depth < 0 | object@GBM$interaction.depth%%1!=0){ cat("\nGBM$interaction.depth must be a positive integer"); test <- FALSE }
           }

           if(!is.numeric(object@GBM$n.minobsinnode)){ cat("\nGBM$n.minobsinnode must be a integer"); test <- FALSE } else{
             if(object@GBM$n.minobsinnode < 0 | object@GBM$n.minobsinnode%%1!=0){ cat("\nGBM$n.minobsinnode must positive "); test <- FALSE }
           }

           if(!is.numeric(object@GBM$shrinkage)){ cat("\nGBM$shrinkage must be a numeric"); test <- FALSE } else{
             if(object@GBM$shrinkage < 0 ){ cat("\nGBM$shrinkage must positive "); test <- FALSE }
           }

           if(!is.numeric(object@GBM$bag.fraction)){ cat("\nGBM$bag.fraction must be a numeric"); test <- FALSE } else{
             if(object@GBM$bag.fraction < 0 | object@GBM$bag.fraction > 1){ cat("\nGBM$bag.fraction must be a 0 to 1 numeric"); test <- FALSE }
           }

           if(!is.numeric(object@GBM$train.fraction)){ cat("\nGBM$train.fraction must be a numeric"); test <- FALSE } else{
             if(object@GBM$train.fraction < 0 | object@GBM$train.fraction > 1){ cat("\nGBM$train.fraction must be a 0 to 1 numeric"); test <- FALSE }
           }

           if(!is.numeric(object@GBM$cv.folds)){ cat("\nGBM$cv.folds must be a integer"); test <- FALSE } else{
             if(object@GBM$cv.folds < 0 | object@GBM$cv.folds%%1!=0){ cat("\nGBM$cv.folds must be a positive integer"); test <- FALSE }
           }

           if(!is.logical(object@GBM$keep.data)){ cat("\nGBM$keep.data must be a logical"); test <- FALSE }

           if(!is.logical(object@GBM$verbose)){ cat("\nGBM$verbose must be a logical"); test <- FALSE }

           #            if(! object@GBM$class.stratify.cv %in% c("bernoulli", "multinomial")){cat("\nGBM$class.stratify.cv must be 'bernoulli' or 'multinomial'") }

           if(! object@GBM$perf.method %in% c('OOB', 'test', 'cv')){cat("\nGBM$perf.method must be 'OOB', 'test' or 'cv'"); test <- FALSE }


           ## GAM ##
           if(! object@GAM$algo %in% c('GAM_mgcv','GAM_gam', 'BAM_mgcv')){cat("\nGAM$algo must be 'GAM_mgcv','GAM_gam' or  'BAM_mgcv'"); test <- FALSE }

           if(! object@GAM$type %in% c('s_smoother','s', 'lo', 'te')){cat("\nGAM$type must be c('s_smoother','s', 'lo' or 'te'"); test <- FALSE }

           if(! is.null(object@GAM$k)){
             if(! is.numeric(object@GAM$k)  ){ cat("\nGAM$k must be a integer"); test <- FALSE } else{
               if(object@GAM$k < -1 | object@GAM$k%%1!=0){ cat("\nGAM$k must be > -1"); test <- FALSE }
             }
           }


           if(!is.numeric(object@GAM$interaction.level)){ cat("\nGAM$interaction.level must be a integer"); test <- FALSE } else{
             if(object@GAM$interaction.level < 0 | object@GAM$interaction.level%%1!=0){ cat("\nGAM$interaction.level must be a positive integer"); test <- FALSE }
           }

           if(!is.null(object@GAM$myFormula)) if(!inherits(object@GAM$myFormula, "formula")){ cat("\nGAM$myFormula must be NULL or a formula object"); test <- FALSE }

           if(!inherits(object@GAM$family, "family")){ cat("\nGAM$family must be a valid family object"); test <- FALSE }

           if(!is.list(object@GAM$control)){cat("\nGAM$control must be a list like that returned by gam.control"); test <- FALSE}
           if(! object@GAM$method %in% c('GCV.Cp','GACV.Cp','REML', 'P-REML', 'ML', 'P-ML')){cat("\nGAM$method must be 'GCV.Cp','GACV.Cp','REML', 'P-REML', 'ML'or 'P-ML'"); test <- FALSE}

           if(sum(! object@GAM$optimizer %in% c('perf','outer', 'newton', 'bfgs', 'optim', 'nlm', 'nlm.fd')) > 0 ){cat("\nGAM$optimizer bad definition (see ?mgcv::gam)") ; test <- FALSE}

           if(!is.logical(object@GAM$select)){ cat("\nGAM$select must be a logical"); test <- FALSE }

           #            knots=NULL,
           #            paraPen=NULL


           ## CTA ##
           if(! object@CTA$method %in% c( 'anova', 'poisson', 'class', 'exp')){cat("\nCTA$method must be 'anova', 'poisson', 'class' or 'exp'"); test <- FALSE }

           #parms = 'default',

           if(!is.list(object@CTA$control)){cat("\nCTA$control must be a list like that returned by rpart.control"); test <- FALSE}
           if(length(object@CTA$cost)){
             if(!is.numeric(object@CTA$cost)){cat("\nCTA$cost must be a non negative cost vector"); test <- FALSE}
           }



           ## ANN ##
           if(!is.numeric(object@ANN$NbCV)){ cat("\nANN$NbCV must be a integer"); test <- FALSE } else{
             if(object@ANN$NbCV < 0 | object@ANN$NbCV%%1!=0){ cat("\nANN$NbCV must be a positive integer"); test <- FALSE }
           }

           if( ( is.null(object@ANN$size) | length(object@ANN$size)>1 ) & object@ANN$NbCV <= 0){ cat("\nANN$size has to be defined as a single integer if ANN$NbCV=0"); test <- FALSE } else{
             if(!is.null(object@ANN$size)) if( !is.numeric(object@ANN$size) | !all( object@ANN$size > 0 ) | !all( object@ANN$size %% 1 == 0 ) ){ cat("\nANN$size must be NULL or a positive (vector of) integer"); test <- FALSE }
           }

           if( ( is.null(object@ANN$decay) | length(object@ANN$decay)>1 ) & object@ANN$NbCV <= 0){ cat("\nANN$decay has to be defined as a single number if ANN$NbCV=0"); test <- FALSE } else{
             if(!is.null(object@ANN$decay)) if( !is.numeric(object@ANN$decay) | !all( object@ANN$decay > 0 ) ){ cat("\nANN$decay must be NULL or a positive (vector of) number"); test <- FALSE }
           }

           if(!is.numeric(object@ANN$rang)){ cat("\nANN$rang must be a numeric"); test <- FALSE } else{
             if(object@ANN$rang < 0 ){ cat("\nANN$rang must be positive"); test <- FALSE }
           }

           if(!is.numeric(object@ANN$maxit)){ cat("\nANN$maxit must be a integer"); test <- FALSE } else{
             if(object@ANN$maxit < 0 | object@ANN$maxit%%1!=0){ cat("\nANN$maxit must be a positive integer"); test <- FALSE }
           }



           ## FDA ##
           if(! object@FDA$method %in% c( 'polyreg', 'mars', 'bruto')){cat("\nFDA$method must be 'polyreg', 'mars' or 'bruto'"); test <- FALSE }
           if(!is.null(object@FDA$add_args)){ if(!is.list(object@FDA$add_args)) {cat("\nFDA$add_args must be a list or NULL"); test <- FALSE } }


           ## SRE ##
           if(!is.numeric(object@SRE$quant)){ cat("\nSRE$quant must be a numeric"); test <- FALSE } else{
             if(object@SRE$quant >= 0.5 | object@SRE$quant < 0){ cat("\nSRE$quant must between 0 and 0.5"); test <- FALSE }
           }



           ## MARS ##
           if(!(object@MARS$type %in% c('simple','quadratic','polynomial','user.defined'))){ cat("\nMARS$type must be 'simple',  'quadratic', 'polynomial' or 'user.defined'"); test <- FALSE}
           if(!is.numeric(object@MARS$interaction.level)){ cat("\nMARS$interaction.level must be a integer"); test <- FALSE } else{
             if(object@MARS$interaction.level < 0 | object@MARS$interaction.level%%1!=0){ cat("\nMARS$interaction.level must be a positive integer"); test <- FALSE }
           }
           if(!is.null(object@MARS$myFormula)) if(!inherits(object@MARS$myFormula, "formula")){ cat("\nMARS$myFormula must be NULL or a formula object"); test <- FALSE }
#            if(!is.numeric(object@MARS$degree)){ cat("\nMARS$degree must be a integer"); test <- FALSE } else{
#              if(object@MARS$degree < 0 | object@MARS$degree%%1!=0){ cat("\nMARS$degree must be a positive integer"); test <- FALSE }
#            }
           if(!is.null(object@MARS$nk)){
             if(object@MARS$nk < 0 | object@MARS$nk%%1!=0){ cat("\nMARS$nk must be a positive integer or NULL if you want to use default parameter"); test <- FALSE }
           }
           if(!is.numeric(object@MARS$penalty)){ cat("\nMARS$penalty must be a integer"); test <- FALSE } else{
             if(object@MARS$penalty < 0 | object@MARS$penalty%%1!=0){ cat("\nMARS$penalty must be a positive integer"); test <- FALSE }
           }
           if(!is.numeric(object@MARS$thresh)){ cat("\nMARS$thresh must be a numeric"); test <- FALSE } else{
             if(object@MARS$thresh < 0 ){ cat("\nMARS$thresh must be positive"); test <- FALSE }
           }
           if(!is.null(object@MARS$nprune)){ if(!is.numeric(object@MARS$nprune)){ cat("\nMARS$nprune must be a numeric or NULL"); test <- FALSE }}
           supported.pmethod <- c('backward', 'none', 'exhaustive', 'forward', 'seqrep', 'cv')
           if(!is.element(object@MARS$pmethod, supported.pmethod)){cat("\nMARS$pmethod must be a one of", supported.pmethod); test <- FALSE }


           ## RF ##
           if(!is.logical(object@RF$do.classif)){ cat("\nRF$do.classif must be a logical"); test <- FALSE }

           if(!is.numeric(object@RF$ntree)){ cat("\nRF$ntree must be a integer"); test <- FALSE } else{
             if(object@RF$ntree < 0 | object@RF$ntree%%1!=0){ cat("\nRF$ntree must be a positive integer"); test <- FALSE }
           }

           if( object@RF$mtry != 'default'){
             if(!is.numeric(object@RF$mtry)){ cat("\nRF$mtry must be a integer"); test <- FALSE } else{
               if(object@RF$mtry < 0 | object@RF$mtry%%1!=0){ cat("\nRF$mtry must be a positive integer"); test <- FALSE }
             }
           }

           if(!is.numeric(object@RF$nodesize)){ cat("\nRF$nodesize must be a integer"); test <- FALSE } else{
             if(object@RF$nodesize < 0 | object@RF$nodesize%%1!=0){ cat("\nRF$nodesize must be a positive integer"); test <- FALSE }
           }

           if(length(object@RF$maxnodes)){
             if(!is.numeric(object@RF$maxnodes)){ cat("\nRF$maxnodes must be a integer"); test <- FALSE } else{
               if(object@RF$maxnodes < 0 | object@RF$maxnodes%%1!=0){ cat("\nRF$maxnodes must be a positive integer"); test <- FALSE }
             }
           }



           ## MAXENT.Phillips ##
           if(!is.character(object@MAXENT.Phillips$path_to_maxent.jar)){ cat("\nMAXENT.Phillips$path_to_maxent.jar must be a character"); test <- FALSE }
           if(!is.null(object@MAXENT.Phillips$memory_allocated)){
             if(!is.numeric(object@MAXENT.Phillips$memory_allocated)){
               cat("\nMAXENT.Phillips$memory_allocated must be a positive integer or NULL for unlimited memory allocation"); test <- FALSE }
           }
           if(!is.character(object@MAXENT.Phillips$background_data_dir)){ cat("\nMAXENT.Phillips$background_data_dir must be 'default' (=> use the same pseudo absences than other models as background) or a path to the directory where your environmental layer are stored"); test <- FALSE }
           tt <- is.character(object@MAXENT.Phillips$maximumbackground) | is.numeric(object@MAXENT.Phillips$maximumbackground)
           if(is.character(object@MAXENT.Phillips$maximumbackground)) if(object@MAXENT.Phillips$maximumbackground != 'default') tt <- FALSE
           if(!tt){ cat("\nMAXENT.Phillips$maximumbackground must be 'default' or numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$maximumiterations)){ cat("\nMAXENT.Phillips$maximumiterations must be a integer"); test <- FALSE } else{
             if(object@MAXENT.Phillips$maximumiterations < 0 | object@MAXENT.Phillips$maximumiterations%%1!=0){ cat("\nMAXENT.Phillips$maximumiterations must be a positive integer"); test <- FALSE }
           }
           if(!is.logical(object@MAXENT.Phillips$visible)){ cat("\nMAXENT.Phillips$visible must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$linear)){ cat("\nMAXENT.Phillips$linear must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$quadratic)){ cat("\nMAXENT.Phillips$quadratic must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$product)){ cat("\nMAXENT.Phillips$product must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$threshold)){ cat("\nMAXENT.Phillips$threshold must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$hinge)){ cat("\nMAXENT.Phillips$hinge must be a logical"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$lq2lqptthreshold)){ cat("\nMAXENT.Phillips$lq2lqptthreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$l2lqthreshold)){ cat("\nMAXENT.Phillips$l2lqthreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$lq2lqptthreshold)){ cat("\nMAXENT.Phillips$lq2lqptthreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$hingethreshold)){ cat("\nMAXENT.Phillips$hingethreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_threshold)){ cat("\nMAXENT.Phillips$beta_threshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_categorical)){ cat("\nMAXENT.Phillips$beta_categorical must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_lqp)){ cat("\nMAXENT.Phillips$beta_lqp must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_hinge)){ cat("\nMAXENT.Phillips$beta_hinge must be a numeric"); test <- FALSE }
		       if(!is.numeric(object@MAXENT.Phillips$betamultiplier)){ cat("\nMAXENT.Phillips$betamultiplier must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$defaultprevalence)){ cat("\nMAXENT.Phillips$defaultprevalence must be a numeric"); test <- FALSE }
           
           ## MAXENT.Phillips.2
           
           ### TO BE DONE ===

#            ## MAXENT.Tsuruoka
# 		       if(!is.numeric(object@MAXENT.Tsuruoka$l1_regularizer)){ cat("\nMAXENT.Tsuruoka$l1_regularizer must be a numeric"); test <- FALSE }
# 		       if(!is.numeric(object@MAXENT.Tsuruoka$l2_regularizer)){ cat("\nMAXENT.Tsuruoka$l2_regularizer must be a numeric"); test <- FALSE }
# 		       if(!is.logical(object@MAXENT.Tsuruoka$use_sgd)){ cat("\nMAXENT.Tsuruoka$use_sgd must be a logical"); test <- FALSE }
# 		       if(!is.numeric(object@MAXENT.Tsuruoka$set_heldout)){ cat("\nMAXENT.Tsuruoka$set_heldout must be a numeric"); test <- FALSE }
# 		       if(!is.logical(object@MAXENT.Tsuruoka$verbose)){ cat("\nMAXENT.Tsuruoka$verbose must be a logical"); test <- FALSE }

           return(test)
         })

##' @rdname BIOMOD.Model.Options-objects
setMethod(
  'show', 
  signature('BIOMOD.Model.Options'),
  function(object){
    .bmCat(" 'BIOMOD.Model.Options' ")
    cat("\n")
    
    ## GLM options
    cat("\nGLM = list( type = '", object@GLM$type, "',", sep="")
    cat("\n            interaction.level = ", object@GLM$interaction.level, ",", sep="")
    cat("\n            myFormula = ",  ifelse(length(object@GLM$myFormula) < 1,'NULL',paste(object@GLM$myFormula[2],object@GLM$myFormula[1],object@GLM$myFormula[3])), ",", sep="")
    cat("\n            test = '", object@GLM$test, "',", sep="")
    cat("\n            family = ", object@GLM$family$family,"(link = '",object@GLM$family$link,"'),", sep="")
    cat("\n            mustart = ", object@GLM$mustart, ",", sep="")
    cat("\n            control = glm.control(", .print.control(object@GLM$control), ") ),", sep="", fill=.Options$width)
    
    ## GBM options
    cat("\n")
    cat("\nGBM = list( distribution = '", object@GBM$distribution, "',", sep="")
    cat("\n            n.trees = ", object@GBM$n.trees, ",", sep="")
    cat("\n            interaction.depth = ", object@GBM$interaction.depth, ",", sep="")
    cat("\n            n.minobsinnode = ", object@GBM$n.minobsinnode, ",", sep="")
    cat("\n            shrinkage = ", object@GBM$shrinkage, ",", sep="")
    cat("\n            bag.fraction = ", object@GBM$bag.fraction, ",", sep="")
    cat("\n            train.fraction = ", object@GBM$train.fraction, ",", sep="")
    cat("\n            cv.folds = ", object@GBM$cv.folds, ",", sep="")
    cat("\n            keep.data = ", object@GBM$keep.data, ",", sep="")
    cat("\n            verbose = ", object@GBM$verbose, ",", sep="")
    #             cat("\n            class.stratify.cv = '", object@GBM$class.stratify.cv, "',", sep="")
    cat("\n            perf.method = '", object@GBM$perf.method, "',", sep="")
    cat("\n            n.cores = ", ifelse(length(object@GBM$n.cores), object@GBM$n.cores,'NULL'), "),", sep="")
    
    ## GAM options
    cat("\n")
    cat("\nGAM = list( algo = '", object@GAM$algo, "',", sep="")
    cat("\n            type = '", object@GAM$type, "',", sep="")
    cat("\n            k = ", ifelse(length(object@GAM$k) < 1,'NULL',object@GAM$k), ",", sep="")
    cat("\n            interaction.level = ", object@GAM$interaction.level, ",", sep="")
    cat("\n            myFormula = ", ifelse(length(object@GAM$myFormula) < 1,'NULL',paste(object@GAM$myFormula[2],object@GAM$myFormula[1],object@GAM$myFormula[3])), ",", sep="")
    cat("\n            family = ", object@GAM$family$family,"(link = '",object@GAM$family$link,"'),", sep="")
    
    if(object@GAM$algo=='GAM_mgcv'){
      cat("\n            method = '", object@GAM$method, "',", sep="")
      cat("\n            optimizer = c('", paste(object@GAM$optimizer,collapse="','"), "'),", sep="")
      cat("\n            select = ", object@GAM$select, ",", sep="")
      cat("\n            knots = ",  ifelse(length(object@GLM$knots) < 1,'NULL',"'user.defined'"), ",", sep="")
      cat("\n            paraPen = ",  ifelse(length(object@GLM$paraPen) < 1,'NULL',"'user.defined'"), ",", sep="")
    }
    
    cat("\n            control = list(", .print.control(object@GAM$control), ") ),", sep="", fill=.Options$width)
    
    
    
    ## CTA options
    cat("\n")
    cat("\nCTA = list( method = '", object@CTA$method, "',", sep="")
    cat("\n            parms = '", object@CTA$parms, "',", sep="")
    cat("\n            cost = ", ifelse(length(object@CTA$cost)<1,'NULL',object@CTA$cost), ",", sep="")
    cat("\n            control = list(", .print.control(object@CTA$control), ") ),", sep="", fill=.Options$width)
    
    ## ANN options
    cat("\n")
    cat("\nANN = list( NbCV = ", object@ANN$NbCV, ",", sep="")
    cat("\n            size = ", ifelse(length(object@ANN$size)<1,'NULL',object@ANN$size), ",", sep="")
    cat("\n            decay = ", ifelse(length(object@ANN$decay)<1,'NULL',object@ANN$decay), ",", sep="")
    cat("\n            rang = ", object@ANN$rang, ",", sep="")
    cat("\n            maxit = ", object@ANN$maxit, "),", sep="")
    
    ## SRE options
    cat("\n")
    cat("\nSRE = list( quant = ", object@SRE$quant, "),", sep="")
    
    ## FDA options
    cat("\n")
    cat("\nFDA = list( method = '", object@FDA$method, "',", sep="")
    cat("\n            add_args = ", ifelse(length(object@FDA$add_args)<1,
                                            'NULL',
                                            paste("list(", paste(.print.control(object@FDA$add_args), collapse=""), ")", sep="")), "),",sep="")
    
    ## MARS options
    cat("\n")
    cat("\nMARS = list( type = '", object@MARS$type, "',", sep="")
    cat("\n             interaction.level = ", object@MARS$interaction.level, ",", sep="")
    cat("\n             myFormula = ",  ifelse(length(object@MARS$myFormula) < 1,'NULL',paste(object@GLM$myFormula[2],object@GLM$myFormula[1],object@GLM$myFormula[3])), ",", sep="")
    #             cat("\n             degree = ", object@MARS$degree, ",", sep="")
    cat("\n             nk = ", ifelse(length(object@MARS$nk) < 1,'NULL',object@MARS$nk), ",", sep="")
    cat("\n             penalty = ", object@MARS$penalty, ",", sep="")
    cat("\n             thresh = ", object@MARS$thresh, ",", sep="")
    cat("\n             nprune = ", ifelse(length(object@MARS$nprune) < 1,'NULL',object@MARS$nprune), ",", sep="")
    cat("\n             pmethod = '", object@MARS$pmethod, "'),", sep="")
    
    ## RF options
    cat("\n")
    cat("\nRF = list( do.classif = ", object@RF$do.classif, ",", sep="")
    cat("\n           ntree = ", object@RF$ntree, ",", sep="")
    cat("\n           mtry = '", object@RF$mtry, "',", sep="")
    cat("\n           nodesize = ", object@RF$nodesize, ",", sep="")
    cat("\n           maxnodes = ", ifelse(length(object@RF$maxnodes) < 1,'NULL',object@RF$maxnodes), "),", sep="")
    
    ## MAXENT.Phillips options
    cat("\n")
    cat("\nMAXENT.Phillips = list( path_to_maxent.jar = '", object@MAXENT.Phillips$path_to_maxent.jar, "',", sep="")
    cat("\n               memory_allocated = ", ifelse(length(object@MAXENT.Phillips$memory_allocated) < 1,'NULL',object@MAXENT.Phillips$memory_allocated), ",", sep="")
    cat("\n               background_data_dir = ", ifelse(is.character(object@MAXENT.Phillips$background_data_dir), "'", ""), object@MAXENT.Phillips$background_data_dir, ifelse(is.character(object@MAXENT.Phillips$background_data_dir), "'", ""), ",", sep="")
    cat("\n               maximumbackground = ", ifelse(is.character(object@MAXENT.Phillips$maximumbackground), "'", ""), object@MAXENT.Phillips$maximumbackground, ifelse(is.character(object@MAXENT.Phillips$maximumbackground), "'", ""), ",", sep="")
    cat("\n               maximumiterations = ", object@MAXENT.Phillips$maximumiterations, ",", sep="")
    cat("\n               visible = ", object@MAXENT.Phillips$visible, ",", sep="")
    cat("\n               linear = ", object@MAXENT.Phillips$linear, ",", sep="")
    cat("\n               quadratic = ", object@MAXENT.Phillips$quadratic, ",", sep="")
    cat("\n               product = ", object@MAXENT.Phillips$product, ",", sep="")
    cat("\n               threshold = ", object@MAXENT.Phillips$threshold, ",", sep="")
    cat("\n               hinge = ", object@MAXENT.Phillips$hinge, ",", sep="")
    cat("\n               lq2lqptthreshold = ", object@MAXENT.Phillips$lq2lqptthreshold, ",", sep="")
    cat("\n               l2lqthreshold = ", object@MAXENT.Phillips$l2lqthreshold, ",", sep="")
    cat("\n               hingethreshold = ", object@MAXENT.Phillips$hingethreshold, ",", sep="")
    cat("\n               beta_threshold = ", object@MAXENT.Phillips$beta_threshold, ",", sep="")
    cat("\n               beta_categorical = ", object@MAXENT.Phillips$beta_categorical, ",", sep="")
    cat("\n               beta_lqp = ", object@MAXENT.Phillips$beta_lqp, ",", sep="")
    cat("\n               beta_hinge = ", object@MAXENT.Phillips$beta_hinge, ",", sep="")
    cat("\n               betamultiplier = ", object@MAXENT.Phillips$betamultiplier, ",", sep="")
    cat("\n               defaultprevalence = ", object@MAXENT.Phillips$defaultprevalence, "),", sep="")
    ## MAXENT.Phillips.2 options
    cat("\n")
    cat("\n MAXENT.Phillips.2 = ")
    cat("\n   list(")
    cat(
      "\n     myFormula = ", 
      cat_formula(object@MAXENT.Phillips.2$myFormula), 
      ",", sep=""
    )
    cat("\n     regmult = ", object@MAXENT.Phillips.2$regmult, ",", sep="")
    cat("\n     regfun = <function>")
    cat("\n   )")
    cat("\n)")
    
    # ## MAXENT.Tsuruoka
    # cat("\n")
    # cat("\nMAXENT.Tsuruoka = list( l1_regularizer = ", object@MAXENT.Tsuruoka$l1_regularizer, ",", sep="")
    # cat("\n                        l2_regularizer = ", object@MAXENT.Tsuruoka$l2_regularizer, ",", sep="")
    # cat("\n                        use_sgd = ", object@MAXENT.Tsuruoka$use_sgd, ",", sep="")
    # cat("\n                        set_heldout = ", object@MAXENT.Tsuruoka$set_heldout, ",", sep="")
    # cat("\n                        verbose = ", object@MAXENT.Tsuruoka$verbose, ")", sep="")
    
    .bmCat()
  })


## TO DO: ===
## moove this function somewhere else
cat_formula <-
  function(
    formula = NULL
  ){
    ifelse(
      length(formula) < 1,
      'NULL',
      paste(
        formula[2],
        formula[1],
        formula[3])
    )
  }


.print.control <- function(ctrl){
  out <-  paste(names(ctrl)[1], " = ", ctrl[[1]], sep="")

  if(length(ctrl) > 1){
    i=2
    while(i <= length(ctrl)){
      if(is.list(ctrl[[i]])){
        out <- c(out, paste(", ", names(ctrl)[i], " = list(",
                            paste(names(ctrl[[i]]), "=",unlist(ctrl[[i]]), sep="", collapse=", "),")", sep=""))
        #         i <- i+length(ctrl[[i]])
        i <- i+1
      } else {
        out <- c(out, paste(", ", names(ctrl)[i], " = ", ctrl[[i]], sep=""))
        i <- i+1
      }
    }
  }
  #   return(toString(out))
  return(out)

}


####################################################################################################
### BIOMOD Storing Results Objects #################################################################
####################################################################################################
setClass("BIOMOD.stored.data",
         representation(inMemory = 'logical',
                        link = 'character'),
         prototype(inMemory=FALSE,
                   link = ''),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.array",
         contains = "BIOMOD.stored.data",
         representation(val = 'array'),
         prototype(val = array()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.data.frame",
         contains = "BIOMOD.stored.data",
         representation(val = 'data.frame'),
         prototype(val = data.frame()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.raster.stack",
         contains = "BIOMOD.stored.data",
         representation(val = 'RasterStack'),
         prototype(val = stack()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.files",
         contains = "BIOMOD.stored.data",
         representation(val = 'character'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.formated.data",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.formated.data'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.models.options",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.Model.Options'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setMethod("load_stored_object", "BIOMOD.stored.data",
          function(obj){
            if(obj@inMemory){
              return(obj@val)
            }

            # different comportement with raster
            if(inherits(obj, "BIOMOD.stored.raster.stack")){
              if( length(obj@link) == 1 & all(grepl(".RData", obj@link)) ){
                return(get(load(obj@link)))
              } else if(all(grepl(".grd", obj@link) | grepl(".img", obj@link))){
                out <- raster::stack(x = obj@link, RAT=FALSE)
                ## rename layer in case of individual projections
                if(all(grepl("individual_projections",obj@link))){
                  # remove directories arch and extention
                  xx <- sub("[:.:].+$", "", sub("^.+individual_projections/", "",obj@link))
                  # remove projection name
                  to_rm <- unique(sub("[^_]+[:_:][^_]+[:_:][^_]+[:_:][^_]+$", "", xx))
                  xx <- sub(to_rm,"", xx)
                  names(out) <- xx
                }
                return(out)
              } # else {
              #                 filesToLoad <- list.files(path=sub("/individual_projections","", obj@link), full.names=T)
              #                 toMatch <- c('.grd$','.img$')
              #                 filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
              #                 if(length(filesToLoad)){
              #                   return(raster::stack(filesToLoad[1]))
              #                 } else {
              #                   filesToLoad <- list.files(path=obj@link, full.names=T)
              #                   toMatch <- c('.grd$','.img$')
              #                   filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
              #                   toMatch <- obj@models.projected
              #                   filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
              #                   proj <- raster::stack(filesToLoad)
              #                   toMatch <- c(obj@link,".img$",'.grd$', .Platform$file.sep)
              #                   names(proj) <- gsub(pattern=paste(toMatch,collapse="|"), "", filesToLoad)
              #                   return(proj)
              #                 }
              #               }
            } else { # for all other stored objects
              return(get(load(obj@link)))
            }
          }

)



setClass(
  "BIOMOD.models.out",
  representation(
    modeling.id = 'character',
    sp.name = 'character',
    expl.var.names = 'character',
    models.computed = 'character',
    models.failed = 'character',
    has.evaluation.data = 'logical',
    rescal.all.models = 'logical',
    models.evaluation = 'BIOMOD.stored.array',
    variables.importances = 'BIOMOD.stored.array',
    models.prediction = 'BIOMOD.stored.array',
    models.prediction.eval = 'BIOMOD.stored.array',
    formated.input.data = 'BIOMOD.stored.formated.data',
    calib.lines = 'BIOMOD.stored.array',
    models.options = 'BIOMOD.stored.models.options',
    link = 'character'
  ),
  prototype(
    modeling.id = as.character(format(Sys.time(), "%s")),
    sp.name = '',
    expl.var.names = '',
    models.computed = '',
    models.failed = '',
    has.evaluation.data = FALSE,
    rescal.all.models = TRUE,
    models.evaluation = new('BIOMOD.stored.array'),
    variables.importances = new('BIOMOD.stored.array'),
    models.prediction = new('BIOMOD.stored.array'),
    models.prediction.eval = new('BIOMOD.stored.array'),
    formated.input.data = new('BIOMOD.stored.formated.data'),
    calib.lines = new('BIOMOD.stored.array'),
    models.options = new('BIOMOD.stored.models.options'),
    link = ''
  ),
  validity = 
    function(object){
      return(TRUE)
    }
)

setClass(
  "BIOMOD.stored.models.out",
  contains = "BIOMOD.stored.data",
  representation(val = 'BIOMOD.models.out'),
  prototype(val = NULL),
  validity = 
    function(object){
      return(TRUE)
    }
)

# #' @rdname BIOMOD.models.out-objects
# #' @docType method
# #' @aliases show, BIOMOD.models.out-method
setMethod(
  'show', 
  signature('BIOMOD.models.out'),
  function(object){
    .bmCat("BIOMOD.models.out")
    cat("\nModeling id :", object@modeling.id, fill=.Options$width)
    cat("\nSpecies modeled :", object@sp.name, fill=.Options$width)
    cat("\nConsidered variables :", object@expl.var.names, fill=.Options$width)
  
    cat("\n\nComputed Models : ", object@models.computed, fill=.Options$width)
    cat("\n\nFailed Models : ", object@models.failed, fill=.Options$width)
    .bmCat()
  }
)


setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

### GETTEURS ###

setMethod("get_predictions", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, evaluation = FALSE){
            # check evaluation data avialability
            if( evaluation & (! obj@has.evaluation.data) ){
              warning("calibration data returned because no evaluation data available")
              evaluation = FALSE
            }

            # select calibration or eval data
            if(evaluation) pred <- obj@models.prediction.eval else pred <- obj@models.prediction

            if(!as.data.frame){
              if(pred@inMemory ){
                return(pred@val)
              } else{
                if(pred@link != ''){
                  return(get(load(pred@link)))
                } else{ return(NULL) }
              }
            } else {
              if(pred@inMemory ){
                mod.pred <- as.data.frame(pred@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE),
                                                 function(x){
                                                   x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                   data.set.id <- x.rev[1]
                                                   cross.valid.id <- x.rev[2]
                                                   algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                   model.id <- paste(obj@sp.name,
                                                                     data.set.id,
                                                                     cross.valid.id,
                                                                     algo.id, sep="_")
                                                   return(model.id)
                                                 }))
                return(mod.pred)
              } else{
                if(pred@link != ''){
                  mod.pred <- as.data.frame(get(load(pred@link)))
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE),
                                                   function(x){
                                                     x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                     data.set.id <- x.rev[1]
                                                     cross.valid.id <- x.rev[2]
                                                     algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                     model.id <- paste(obj@sp.name,
                                                                       data.set.id,
                                                                       cross.valid.id,
                                                                       algo.id, sep="_")
                                                     return(model.id)
                                                   }))
                  return(mod.pred)
                } else{ return(NULL) }
              }

            }
          }
)

setMethod("get_evaluations", "BIOMOD.models.out",
          function(obj, ...){
            args <- list(...)

            ## fill some additional parameters
            as.data.frame <- ifelse(!is.null(args$as.data.frame), args$as.data.frame, FALSE)

            out <- NULL
            if(obj@models.evaluation@inMemory ){
              out <- obj@models.evaluation@val
            } else{
              if(obj@models.evaluation@link != ''){
                out <- get(load(obj@models.evaluation@link))
              }
            }

            ## transform into data.frame object if needed
            if(as.data.frame){
              tmp <- reshape::melt.array(out,varnames=c("eval.metric", "test","m","r","d"))
              model_names <- unique(apply(tmp[,c("m","r","d"), drop=F], 1, paste, collapse="_"))
              out <- data.frame() #NULL
              for(mod in model_names){
                m = unlist(strsplit(mod,"_"))[1]
                r = unlist(strsplit(mod,"_"))[2]
                d = unlist(strsplit(mod,"_"))[3]
                eval.met = as.character(unique(tmp[which( tmp$m == m & tmp$r == r & tmp$d == d), "eval.metric", drop=T]))
                for(em in eval.met){
                  out <- rbind(out,
                               data.frame( Model.name = mod,
                                           Eval.metric = em,
                                           Testing.data = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Testing.data"), "value", drop=T]),
                                           Evaluating.data = ifelse("Evaluating.data" %in% tmp$test, as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Evaluating.data"), "value", drop=T]), NA ),
                                           Cutoff = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Cutoff"), "value", drop=T]),
                                           Sensitivity = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Sensitivity"), "value", drop=T]),
                                           Specificity = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Specificity"), "value", drop=T]) )
                  )
                } # end loop on eval metric
              } # end loop on models names

            }

            return(out)
          }
)


setMethod("get_calib_lines", "BIOMOD.models.out",
   function(obj, as.data.frame = FALSE, ...){
     calib_lines <- load_stored_object(obj@calib.lines)
     return(calib_lines)
   }

)


setMethod("get_variables_importance", "BIOMOD.models.out",
          function(obj, ...){
            if(obj@variables.importances@inMemory ){
              return(obj@variables.importances@val)
            } else{
              if(obj@variables.importances@link != ''){
                return(get(load(obj@variables.importances@link)))
              } else{ return(NA) }
            }
          }
)



setMethod("get_options", "BIOMOD.models.out",
          function(obj){
            if(obj@models.options@inMemory ){
              return(obj@models.options@val)
            } else{
              if(obj@models.options@link != ''){
                return(get(load(obj@models.options@link)))
              } else{ return(NA) }
            }
          }
)

setMethod("get_formal_data", "BIOMOD.models.out",
          function(obj, subinfo = NULL){
            if(is.null(subinfo)){
              if(obj@formated.input.data@inMemory ){
                return(obj@formated.input.data@val)
              } else{
                if(obj@formated.input.data@link != ''){
                  data <- get(load(obj@formated.input.data@link))
                  return(data)
                } else{ return(NA) }
              }
            } else if(subinfo == 'MinMax'){
              return(apply(get_formal_data(obj, "expl.var"),2, function(x){
                if(is.numeric(x)){
                  return( list(min = min(x,na.rm=T), max = max(x, na.rm=T) ) )
                } else if(is.factor(x)){
                  return(list(levels = levels(x)))
                }
              }) )
            } else if(subinfo == 'expl.var'){
              return(as.data.frame(get_formal_data(obj)@data.env.var))
            } else if(subinfo == 'expl.var.names'){
              return(obj@expl.var.names)
            } else if(subinfo == 'resp.var'){
              return(as.numeric(get_formal_data(obj)@data.species))
            } else if(subinfo == 'eval.resp.var'){
              return(as.numeric(get_formal_data(obj)@eval.data.species))
            } else if(subinfo == 'eval.expl.var'){
              return(as.data.frame(get_formal_data(obj)@eval.data.env.var))
            } else{
              stop("Unknown subinfo tag")
            }

          }
)

setMethod("get_built_models", "BIOMOD.models.out",
          function(obj, ...){
            return(obj@models.computed)
          }
)




####################################################################################################
### BIOMOD Storing Projection Objects ##############################################################
####################################################################################################
# setClass("BIOMOD.projection",
#          representation(proj.names = 'character',
#                         sp.name = 'character',
#                         expl.var.names = 'character',
#                         models.computed = 'character',
#                         models.failed = 'character',
#                         models.thresholds = 'BIOMOD.stored.array',
#                         models.prediction = 'BIOMOD.stored.array',
#                         formated.input.data = 'BIOMOD.stored.formated.data',
#                         calib.lines = 'BIOMOD.stored.array',
#                         models.options = 'BIOMOD.stored.models.options'),
#          prototype(sp.name='',
#                    expl.var.names = '',
#                    models.computed='',
#                    models.failed='',
#                    models.evaluation = new('BIOMOD.stored.array'),
#                    variables.importances = new('BIOMOD.stored.array'),
#                    models.prediction = new('BIOMOD.stored.array'),
#                    formated.input.data = new('BIOMOD.stored.formated.data'),
#                    calib.lines = new('BIOMOD.stored.array'),
#                    models.options = new('BIOMOD.stored.models.options')),
#          validity = function(object){
#            return(TRUE)
#            })

setClass("BIOMOD.projection.out",
         representation(proj.names = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.projected = 'character',
                        scaled.models = 'logical',
                        modeling.object = 'BIOMOD.stored.data',
                        modeling.object.id = 'character',
                        type = 'character',
                        proj = 'BIOMOD.stored.data',
                        xy.coord = 'matrix'),
         prototype(proj.names = '',
                   sp.name='',
                   expl.var.names='',
                   models.projected='',
                   scaled.models=TRUE,
                   modeling.object.id='',
                   type='',
                   xy.coord = matrix()),
         validity = function(object){
           return(TRUE)
         })

setMethod("get_projected_models", "BIOMOD.projection.out",
          function(obj){
            return(obj@models.projected)
          })

setMethod("get_predictions", "BIOMOD.projection.out",
          function(obj, as.data.frame=FALSE, full.name=NULL, model=NULL, run.eval=NULL, data.set=NULL){
            models_selected <- get_projected_models(obj)
            if(length(full.name)){
              models_selected <- intersect(full.name, models_selected)
            } else if(length(model) | length(run.eval) | length(data.set)){
              # models subselection according to model, run.eval and sata.set parameters
              if(length(model)) grep_model <- paste("(",paste(model,collapse="|"),")", sep="") else grep_model = "*"
              if(length(run.eval)) grep_run.eval <- paste("(",paste(run.eval,collapse="|"),")", sep="") else grep_run.eval = "*"
              if(length(data.set)) grep_data.set <- paste("(",paste(data.set,collapse="|"),")", sep="") else grep_data.set = "*"
              grep_full <- paste(grep_data.set,"_",grep_run.eval,"_",grep_model,"$",sep="")

              models_selected <- grep(pattern=grep_full, models_selected, value=T)
            }

            #             cat("\n*** models_selected = ", models_selected)

            if (length(models_selected)){
              proj <- load_stored_object(obj@proj)
              names(proj) <- obj@models.projected
              if(inherits(proj,'Raster')){
                proj <- raster::subset(proj,models_selected,drop=FALSE)
              } else {
                if(length(dim(proj)) == 4){ ## 4D arrays
                  proj <- proj[,.extractModelNamesInfo(model.names=models_selected,info='models'),
                               .extractModelNamesInfo(model.names=models_selected,info='run.eval'),
                               .extractModelNamesInfo(model.names=models_selected,info='data.set'), drop=FALSE]
                } else{ ## matrix (e.g. from ensemble models projections)
                  proj <- proj[,models_selected,drop=FALSE]
                }

              }

              if(as.data.frame){
                proj <- as.data.frame(proj)
                ## set correct names

                if(obj@type == 'array' & sum(!(names(proj) %in% models_selected))>0 ){ ## from array & not valid names
                  names(proj) <- unlist(lapply(strsplit(names(proj),".", fixed=TRUE),
                                                   function(x){
                                                     x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                     data.set.id <- x.rev[1]
                                                     cross.valid.id <- x.rev[2]
                                                     algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                     model.id <- paste(obj@sp.name,
                                                                       data.set.id,
                                                                       cross.valid.id,
                                                                       algo.id, sep="_")
                                                     return(model.id)
                                                   }))
                }

                # reorder the data.frame
                proj <- proj[,models_selected]

              }
            } else{
              proj <- NULL
            }

            return(proj)
          })




# 2.3 other functions
setMethod('plot', signature(x='BIOMOD.projection.out', y="missing"),
          function(x,col=NULL, str.grep=NULL){
            models_selected <- x@models.projected
            if(length(str.grep)){
              models_selected <- grep(paste(str.grep,collapse="|"), models_selected,value=T)
            }

            if(!length(models_selected)) stop("invalid str.grep arg")



            if(inherits(x@proj, "BIOMOD.stored.raster.stack")){
              requireNamespace("rasterVis")

              ## define the breaks of the color key
              my.at <- seq(0,1000,by=100)
              ## the labels will be placed vertically centered
              my.labs.at <- seq(0,1000,by=250)
              ## define the labels
              my.lab <- seq(0,1000,by=250)
              ## define colors
              #               my.col <- colorRampPalette(c("red4","orange4","yellow4","green4"))(100)
              my.col <- colorRampPalette(c("grey90","yellow4","green4"))(100)

              ## try to use levelplot function
              try_plot <- try(
                levelplot(get_predictions(x, full.name=models_selected),
                          at=my.at, margin=T, col.regions=my.col,
                          main=paste(x@sp.name,x@proj.names,"projections"),
                          colorkey=list(labels=list(
                            labels=my.lab,
                            at=my.labs.at)))
              )
              if(! inherits(try_plot,"try-error")){ ## produce plot
                print(try_plot)
              } else{## try classical plot
                cat("\nrasterVis' levelplot() function failed. Try to call standard raster plotting function.",
                    "It can lead to unooptimal representations.",
                    "You should try to do it by yourself extracting predicions (see : get_predictions() function)", fill=options()$width)
                try_plot <- try(
                  plot(get_predictions(x, full.name=models_selected))
                )
              }

              if(inherits(try_plot,"try-error")){ # try classical plot
                cat("\n Plotting function failed.. You should try to do it by yourself!")
              }

            } else if(inherits(x@proj, "BIOMOD.stored.array")){
              if(ncol(x@xy.coord) != 2){
                cat("\n ! Impossible to plot projections because xy coordinates are not available !")
              } else {
                multiple.plot(Data = get_predictions(x, full.name=models_selected, as.data.frame=T), coor = x@xy.coord)
              }

            } else {cat("\n !  Biomod Projection plotting issue !", fill=.Options$width)}

          })

##' @rdname BIOMOD.projection.out-objects
##' @docType method
##' @aliases show, BIOMOD.projection.out-method
setMethod('show', signature('BIOMOD.projection.out'),
          function(object){
            .bmCat("'BIOMOD.projection.out'")
            cat("\nProjection directory :", paste(object@sp.name,"/",object@proj.names, sep=""), fill=.Options$width)
            cat("\n")
            cat("\nsp.name :", object@sp.name, fill=.Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill=.Options$width)
            cat("\n")
            cat("\nmodeling id :", object@modeling.object.id ,"(",object@modeling.object@link ,")", fill=.Options$width)
            cat("\nmodels projected :", toString(object@models.projected), fill=.Options$width)

            .bmCat()
          })


setMethod(f='free',
          signature='BIOMOD.projection.out',
          definition = function(obj){
            if(inherits(obj@proj,"BIOMOD.stored.array")){
              obj@proj@val = array()
            } else if(inherits(obj@proj,"BIOMOD.stored.raster.stack")){
              obj@proj@val = stack()
            } else{
              obj@proj@val = NULL
            }
            obj@proj@inMemory = FALSE

            return(obj)
          })

####################################################################################################
### BIOMOD Storing Ensemble Modeling Objects #######################################################
####################################################################################################

##' @name BIOMOD.EnsembleModeling.out-class
##' @rdname BIOMOD.EnsembleModeling.out-objects
##' @docType class
##' 
##' @aliases BIOMOD.EnsembleModeling.out-class
##' @aliases BIOMOD.EnsembleModeling.out
##' 
##' @title BIOMOD_EnsembleModeling() outputs objects class
##' 
##' @description
##' EnsembleModeling objects are created, used and returned by BIOMOD
##' functions. It's contains information relative to an \pkg{biomod2}
##' ensemble modeling procedure.
##' 
##' - output of: \code{\link[biomod2]{BIOMOD_EnsembleModeling}}
##' - input of: \code{\link[biomod2]{BIOMOD_EnsembleForecasting}}
##' 
##' @slot sp.name "character", species name
##' @slot expl.var.names "character", explanatory variables
##' names
##' @slot models.out.obj "BIOMOD.stored.models.out", object which
##' contains information on individuals models that have been combined
##' @slot eval.metric "character", evaluation metrics chose for
##' models selection
##' @slot eval.metric.quality.threshold "numeric", thresholds
##' defined for models selection 
##' @slot em.computed "character", ensemble models built names
##' @slot em.by "character", way models are combined
##' @slot em.models "ANY", list of built biomod2.ensemble.models
##' objects
##' @slot modeling.id "character", the id of the whole
##' modelling process
##' @slot link "character", the path to corresponding hard drive
##' saved object
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_Projection}}, 
##' \code{\link[biomod2]{BIOMOD_Modeling}}, 
##' \code{\link[biomod2]{BIOMOD_EnsembleModeling}}, 
##' \code{\link[biomod2]{BIOMOD_EnsembleForecasting}}
##' 
##' @keywords models
##' @keywords ensemble
##' @author Damien Georges 
##' 
##' @examples
##' showClass("BIOMOD.EnsembleModeling.out")
##' 
setClass("BIOMOD.EnsembleModeling.out",
         representation(sp.name = 'character',
                        expl.var.names = 'character',
                        models.out.obj = 'BIOMOD.stored.models.out',
                        eval.metric = 'character',
                        eval.metric.quality.threshold = 'numeric',
                        em.computed = 'character',
                        em.by = 'character',
                        em.models = 'ANY',
                        modeling.id = 'character',
                        link = 'character'),
         #                         em.models.kept = 'list',
         #                         em.prediction = 'BIOMOD.stored.array',
         #                         em.evaluation = 'BIOMOD.stored.array',
         #                         em.res = 'list',
         #                         em.ci.alpha = 'numeric',
         #                         em.weight = 'list',
         #                         em.bin.tresh = 'list'),
         prototype( sp.name = '',
                    expl.var.names = '',
                    models.out.obj = new('BIOMOD.stored.models.out'),
                    eval.metric = '',
                    eval.metric.quality.threshold = 0,
                    em.models = NULL,
                    em.computed = character(),
                    modeling.id = '.',
                    link = ''),
         #                     em.models.kept = NULL,
         #                     em.prediction = NULL,
         #                     #                     em.evaluation = NULL,
         #                     em.res = list(),
         #                     em.ci.alpha = 0.05,
         #                     em.weight = list(),
         #                     em.bin.tresh = list()),
         validity = function(object){
           return(TRUE)
         })

##' @rdname BIOMOD.EnsembleModeling.out-objects
##' @docType method
##' @param object a BIOMOD.EnsembleModeling.out object
setMethod('show', signature('BIOMOD.EnsembleModeling.out'),
          function(object){
            .bmCat("'BIOMOD.EnsembleModeling.out'")
            cat("\nsp.name :", object@sp.name, fill=.Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill=.Options$width)
            cat("\n")
            cat("\nmodels computed:", toString(object@em.computed), fill=.Options$width)

            .bmCat()
          })

setMethod("get_needed_models", "BIOMOD.EnsembleModeling.out",
          function(obj, subset='all', ...){
            add.args <- list(...)
            needed_models <- lapply(obj@em.models, function(x){
              return(x@model)
            })
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          }
)

setMethod("get_needed_models", "BIOMOD.EnsembleModeling.out",
          function(obj, selected.models='all', ...){           ### MODIFIED ROBIN : ici j'ai renom? "subset" en "selected.models" pour que ?a soit le m?me nom que dans les autres fonctions
            add.args <- list(...)
            if(selected.models[[1]]=="all") selected.index <- c(1:length(obj@em.models)) else selected.index <- which(names(obj@em.models) %in% selected.models) ### MODIFIED ROBIN : ici je selectionne uniquement le sous-ensemble des mod?les que l'utilisateur a sp?cifi?.
            needed_models <- lapply(obj@em.models[selected.index], function(x) return(x@model))  ### MODIFIED ROBIN : ici je selectionne uniquement le sous-ensemble des mod?les que l'utilisateur a sp?cifi?.
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          }
)


setMethod("get_kept_models", "BIOMOD.EnsembleModeling.out",
          function(obj, model, ...){
            if(is.character(model) | is.numeric(model)){
              return(obj@em.models[[model]]@model)
            } else{
              kept_mod <- lapply(obj@em.models, function(x){return(x@model)})
              names(kept_mod) <- names(obj@em.models)
              return(kept_mod)
            }

          }
)

setMethod("get_evaluations", "BIOMOD.EnsembleModeling.out",
          function(obj, ...){
            args <- list(...)

            ## fill some additional parameters
            as.data.frame <- ifelse(!is.null(args$as.data.frame), args$as.data.frame, FALSE)

            out <- list()

            ## list of computed models
            models <- obj@em.computed

            ## extract evaluation scores as a list
            for(mod in models){
              out[[mod]] <- obj@em.models[[mod]]@model_evaluation[,,drop=F]
            }

            ## transform into data.frame object if needed
            if(as.data.frame){
              tmp <- melt(out, varnames=c("eval.metric", "test"))
              tmp$model.name <- sapply(tmp$L1, function(x){paste(unlist(strsplit(x, "_"))[-1], collapse="_")})
              out <- data.frame() #NULL
              for(mod in unique(tmp$model.name)){
                eval.met = as.character(unique(tmp[which( tmp$model.name == mod ), "eval.metric", drop=T]))
                for(em in eval.met){
                  out <- rbind(out,
                               data.frame( Model.name = mod,
                                           Eval.metric = em,
                                           Testing.data = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Testing.data"), "value", drop=T]),
                                           Evaluating.data = ifelse("Evaluating.data" %in% tmp$test, as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Evaluating.data"), "value", drop=T]), NA ),
                                           Cutoff = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Cutoff"), "value", drop=T]),
                                           Sensitivity = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Sensitivity"), "value", drop=T]),
                                           Specificity = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Specificity"), "value", drop=T]) )
                  )
                } # end loop on eval metric
              } # end loop on models names
            } # end as.data.frame == TRUE

            return(out)
          }
)


setMethod("get_variables_importance", "BIOMOD.EnsembleModeling.out",
          function(obj, ...){
            vi <- NULL
            for (mod in get_built_models(obj) ){
              (vi_tmp <- obj@em.models[[mod]]@model_variables_importance)
              vi <- abind::abind(vi, vi_tmp, along=3)
            }
            dimnames(vi)[[3]] <- get_built_models(obj)

            return(vi)
          }
)

setMethod("get_built_models", "BIOMOD.EnsembleModeling.out",
          function(obj, ...){
            return(obj@em.computed)
          })

setMethod("get_predictions", "BIOMOD.EnsembleModeling.out",
  function(obj, ...){
    ## note: ensemble models predicitons are stored within the directory
    ##  <sp.name>/.BIOMOD_DATA/<modelling.id>/ensemble.models/ensemble.models.projections/
    ##  This function is just a friendly way to load this data

    ## get the path to projections files we want to load
    files.to.load <- file.path(obj@sp.name, ".BIOMOD_DATA", obj@modeling.id, "ensemble.models",
                              "ensemble.models.predictions", paste0(obj@em.computed, ".predictions"))
    ## load and merge projection files within a data.frame
    bm.pred <- do.call(cbind, lapply(files.to.load, function(ftl) get(load(ftl))))
    colnames(bm.pred) <- obj@em.computed
    return(bm.pred)
  }
)

"/home/georgeda/Work/BIOMOD/RForge/tests/workdir/GuloGulo/.BIOMOD_DATA/test//GuloGulo_EMmeanByTSS_mergedAlgo_mergedRun_mergedData.predictions"
####################################################################################################
### BIOMOD Storing Ensemble Forecasting Objects ####################################################
####################################################################################################

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# if( !isGeneric( ".Models.prepare.data" ) ) {
# }

setMethod(
  '.Models.prepare.data', 
  signature(
    data = 'BIOMOD.formated.data'
  ),
  function(
    data, 
    NbRunEval, 
    DataSplit, 
    Yweights = NULL, 
    Prevalence = NULL, 
    do.full.models = TRUE, 
    DataSplitTable = NULL
  ){
    list.out <- list()
    name <- paste0(data@sp.name, '_AllData')
    xy <- data@coord
    # dataBM <- data.frame(cbind(data@data.species, data@data.env.var))
    # colnames(dataBM)[1] <- data@sp.name
    dataBM <- 
      bind_cols(
        tibble(
          !!data@sp.name := data@data.species
        ),
        data@data.env.var
      )
    
    # dealing with evaluation data
    if(data@has.data.eval){
      evalDataBM <- data.frame(cbind(data@eval.data.species,data@eval.data.env.var[,,drop=FALSE]))
      colnames(evalDataBM)[1] <- data@sp.name
      eval.xy <- data@eval.coord
    } else{ evalDataBM <- eval.xy <- NULL }

    ### Calib/Valid lines
    if(!is.null(DataSplitTable)){
      calibLines <- DataSplitTable
      colnames(calibLines) <- paste('_RUN',1:ncol(calibLines), sep='')
    } else {
      if(NbRunEval == 0){ # take all available data
        calibLines <- matrix(rep(TRUE,length(data@data.species)),ncol=1)
        colnames(calibLines) <- '_Full'
      } else {
        calibLines <- .SampleMat(data.sp = data@data.species,
                                 dataSplit = DataSplit,
                                 nbRun = NbRunEval,
                                 data.env = data@data.env.var)
        if(do.full.models){
          calibLines <- cbind(calibLines, rep(TRUE,length(data@data.species)))
          colnames(calibLines)[NbRunEval+1] <- '_Full'
        }
      }
    }
    ## force calib.lines object to be 3D array
    if(length(dim(calibLines)) < 3 ){
      dn_tmp <- dimnames(calibLines) ## keep track of dimnames
      dim(calibLines) <- c(dim(calibLines),1)
      dimnames(calibLines) <- list(dn_tmp[[1]], dn_tmp[[2]], "_AllData")
    }

    if(is.null(Yweights)){ # 1 for all points
      if(!is.null(Prevalence)){
        cat("\n\t> Automatic weights creation to rise a", Prevalence,"prevalence")
        Yweights <- .automatic_weights_creation(data@data.species ,prev=Prevalence)
      } else{
        cat("\n\t> No weights : all observations will have the same weight")
        Yweights <- rep(1,length(data@data.species))
      }

    }
    list.out[[name]] <- list(name=name,
                             xy=xy,
                             dataBM=dataBM,
                             calibLines=calibLines,
                             Yweights = Yweights,
                             evalDataBM = evalDataBM,
                             eval.xy = eval.xy)
    return(list.out)
  })

setMethod('.Models.prepare.data', signature(data='BIOMOD.formated.data.PA'),
          function(data, NbRunEval, DataSplit, Yweights=NULL, Prevalence=NULL, do.full.models=TRUE, DataSplitTable=NULL){
            list.out <- list()
            formal_weights <- Yweights
            for(pa in 1:ncol(data@PA)){
              Yweights <- formal_weights
              name <- paste(data@sp.name,"_",colnames(data@PA)[pa],sep="")
              xy <- data@coord[data@PA[,pa],]
              resp <- data@data.species[data@PA[,pa]] # response variable (with pseudo absences selected)
              resp[is.na(resp)] <- 0
              dataBM <- data.frame(cbind(resp,
                                         data@data.env.var[data@PA[,pa],,drop=FALSE]))
              colnames(dataBM)[1] <- data@sp.name

              ### Calib/Valid lines
              if(!is.null(DataSplitTable)){
                if(length(dim(DataSplitTable))==2){
                  calibLines <- DataSplitTable
                } else {
                  calibLines <- asub(DataSplitTable,pa,3,drop=TRUE)
                }
                colnames(calibLines) <- paste('_RUN',1:ncol(calibLines), sep='')
                calibLines[which(!data@PA[,pa]),] <- NA
              } else{
                if(NbRunEval == 0){ # take all available data
                  calibLines <- matrix(NA,nrow=length(data@data.species),ncol=1)
                  calibLines[data@PA[,pa],1] <- TRUE
                  colnames(calibLines) <- '_Full'
                } else {
                  calibLines <- matrix(NA,nrow=length(data@data.species),ncol=NbRunEval)
                  sampled.mat <- .SampleMat(data.sp = data@data.species[data@PA[,pa]],
                                            dataSplit = DataSplit,
                                            nbRun = NbRunEval,
                                            data.env = data@data.env.var[data@PA[,pa], , drop = FALSE])
                  calibLines[data@PA[,pa],] <- sampled.mat
                  colnames(calibLines) <- colnames(sampled.mat)
                  if(do.full.models){
                    calibLines <- cbind(calibLines, rep(NA,length(data@data.species)))
                    calibLines[data@PA[,pa],NbRunEval+1] <- TRUE
                    colnames(calibLines)[NbRunEval+1] <- '_Full'
                  }
                }
              }

              ## force calib.lines object to be 3D array
              if(length(dim(calibLines)) < 3 ){
                dn_tmp <- dimnames(calibLines) ## keep track of dimnames
                dim(calibLines) <- c(dim(calibLines),1)
                dimnames(calibLines) <- list(dn_tmp[[1]], dn_tmp[[2]], paste("_PA",pa, sep=""))
              }

              # dealing with evaluation data
              if(data@has.data.eval){
                evalDataBM <- data.frame(cbind(data@eval.data.species,data@eval.data.env.var))
                colnames(evalDataBM)[1] <- data@sp.name
                eval.xy <- data@eval.coord
              } else{ evalDataBM <- eval.xy <- NULL }

              if(is.null(Yweights)){ # prevalence of 0.5... may be parametrize
                if(is.null(Prevalence)) Prevalence <- 0.5

                cat("\n\t\t\t! Weights where automatically defined for", name, "to rise a", Prevalence, "prevalence !")


                Yweights <- rep(NA, length(data@data.species))
                Yweights[data@PA[,pa]] <- .automatic_weights_creation(as.numeric(dataBM[,1]) ,prev=Prevalence)#, subset=data@PA[,pa])
              } else{
                # remove useless weights
                Yweights[!data@PA[,pa]] <- NA
              }

              list.out[[name]] <- list(name=name,
                                       xy=xy,
                                       dataBM=dataBM,
                                       calibLines=calibLines,
                                       Yweights = Yweights,
                                       evalDataBM = evalDataBM,
                                       eval.xy = eval.xy)
            }
            return(list.out)
          })


.automatic_weights_creation <- function(resp,prev=0.5, subset=NULL){
  if(is.null(subset)) subset<- rep(TRUE, length(resp))

  nbPres <- sum(resp[subset], na.rm=TRUE)
  nbAbsKept <- sum(subset, na.rm=T) - sum(resp[subset], na.rm=TRUE) # The number of true absences + pseudo absences to maintain true value of prevalence
  Yweights <- rep(1,length(resp))

  if(nbAbsKept > nbPres){ # code absences as 1
    Yweights[which(resp>0)] <- (prev * nbAbsKept) / (nbPres * (1-prev))
  } else{ # code presences as 1
    Yweights[which(resp==0 | is.na(resp))] <- (nbPres * (1-prev)) / (prev * nbAbsKept)
  }
  Yweights = round(Yweights[])
  Yweights[!subset] <- 0

  return(Yweights)
}


# ##' @name BIOMOD.models.out-RemoveProperly
# ##' @aliases RemoveProperly
# ##' @aliases RemoveProperly, BIOMOD.models.out-method
# ##' @title remove properly BIOMOD_Modeling outputs
# ##' @description
# ##' Functions to free properly a \code{\link[biomod2]{BIOMOD_Modeling}}
# ##' outputs
# ##' 
# ##' @param obj \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} 
# ##'   object
# ##' @param obj.name the name of object in current environment, 
# ##'   automatically filled
# ##' @param ... extra arguments (not implemented yet)
# ##' 
# ##' @details
# ##' This function will remove all objects created during a \code{biomod2}
# ##' modeling run. It will free both objects saved in memory and objects
# ##' saved on hard drive.
# ##' @author Wilfried Thuiller, Damien Georges
# ##' @examples
# ##'   \dontrun{
# ##'     ##' species occurrences
# ##'     DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
# ##'                                         package="biomod2"), row.names = 1)
# ##'     head(DataSpecies)
# ##'     
# ##'     ##' the name of studied species
# ##'     myRespName <- 'VulpesVulpes'
# ##'     
# ##'     ##' the presence/absences data for our species 
# ##'     myResp <- as.numeric(DataSpecies[,myRespName])
# ##'     
# ##'     ##' the XY coordinates of species data
# ##'     myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
# ##'     
# ##'     
# ##'     ##' Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
# ##'     myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd", 
# ##'                                          package="biomod2"),
# ##'                             system.file( "external/bioclim/current/bio4.grd", 
# ##'                                          package="biomod2"), 
# ##'                             system.file( "external/bioclim/current/bio7.grd", 
# ##'                                          package="biomod2"),  
# ##'                             system.file( "external/bioclim/current/bio11.grd", 
# ##'                                          package="biomod2"), 
# ##'                             system.file( "external/bioclim/current/bio12.grd", 
# ##'                                          package="biomod2"))
# ##'     
# ##'     ##' 1. Formatting Data
# ##'     myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
# ##'                                          expl.var = myExpl,
# ##'                                          resp.xy = myRespXY,
# ##'                                          resp.name = myRespName)
# ##'     
# ##'     ##' 2. Defining Models Options using default options.
# ##'     myBiomodOption <- BIOMOD_ModelingOptions()
# ##'     
# ##'     ##' 3. Doing Modelisation
# ##'     
# ##'     myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
# ##'                                          models = c('SRE'), 
# ##'                                          models.options = myBiomodOption, 
# ##'                                          NbRunEval=1, 
# ##'                                          DataSplit=80, 
# ##'                                          Prevalence=0.5, 
# ##'                                          VarImport=0, 
# ##'                                          models.eval.meth = c('TSS','ROC'),
# ##'                                          do.full.models=FALSE,
# ##'                                          modeling.id="test2")
# ##'     
# ##'     ##' files have been created on hard drive
# ##'     list.files(myRespName,all.files=TRUE,recursive=TRUE)
# ##'     
# ##'     ##' remove properly the modeling objects and all the file saved on hard drive
# ##'     RemoveProperly(myBiomodModelOut)
# ##'     
# ##'     ##' check files had been removed
# ##'     list.files(myRespName,all.files=TRUE,recursive=TRUE)
# ##'   }
# ##'   
# ##' @export
# ##' @docType methods
# ##' @rdname RemoveProperly-methods
# setMethod(
#   "RemoveProperly", "BIOMOD.models.out",
#   function(obj, obj.name=deparse(substitute(obj)))
#   {
#     cat("\n\t> Removing .BIOMOD_DATA files...")
#     unlink(file.path(obj@sp.name, ".BIOMOD_DATA", obj@modeling.id), recursive=T, force=TRUE)
#     cat("\n\t> Removing models...")
#     unlink(file.path(obj@sp.name, "models", obj@modeling.id), recursive=T, force=TRUE)
#     cat("\n\t> Removing object hard drive copy...")
#     unlink(obj@link, recursive=T, force=TRUE)
#     cat("\n\t> Removing object from memory...")
#     rm(list=obj.name,envir=sys.frame(-2))
#     cat("\nCompleted!")
#   }
# )
##' @name BIOMOD_ConvertOldRun
##' @title Convert objects and outputs from BIOMOD.xx into biomod2.xx
##' objects and outputs
##' @description 
##' This function converts workspace, modelling outputs, results and
##' objects created with version xx of \pkg{BIOMOD} into \pkg{biomod2}
##' objects and re-organized the directories to be used with \pkg{biomod2}
##' @param savedObj a BIOMOD.1.xx workspace image. It's a .Rdata file
##'   named 'Biomod_run.RData' for plurispecific run or 
##'   Your_Species_Name.RData if you have done monospecific modelling
##' @param path Optional path to your 'savedObj' if you don't have give
##'   the full path to the object
##'   
##' @details 
##' This function is useful to convert former \pkg{BIOMOD} runs into the
##' new \pkg{biomod2} object structure. This is mostly interesting in the
##' case users want to relaunched some projections or analyses within the
##' \pkg{biomod2} new structure. Returned 'BIOMOD.models.out' objects can
##' be then used as classic object for making projections for instance
##' (\code{\link{BIOMOD_Projection}}).
##' 
##' Be aware that because \pkg{biomod2} has strongly changed between the
##' first and second version, some new additional functions and 
##' information could not be used with converted objects (i.e. Calibration
##' Lines access, Maxent run, SRE projections...).   
##' 
##' @return A list of 'BIOMOD.models.out' (one per species modeled)
##' containing information of your old run. Specific directories are also
##' created on your hard drive (see \code{\link{BIOMOD_Modeling}}) 
##' 
##' @author Damien Georges
##' @seealso \code{\link{BIOMOD_Modeling}}
##' @keywords models
##' 
BIOMOD_ConvertOldRun <- function(savedObj, path = NULL){
  .bmCat("BIOMOD results migration")

  compress.arg = TRUE #ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  # 1. Check path exists and all objects needed exists too -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(!file.exists(savedObj)){
    stop("Input object doesn't exist")
  }

  if(is.null(path)){
    path = sub(tail(unlist(strsplit(savedObj,'/')),1), '', savedObj)
  } else{ # add / at the end of path
    if(tail(unlist(strsplit(path)),1) != '/'){
      path = paste(path,"/",sep="")
    }
  }

  # 2. Keep image of current workSpace -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   save.image("WS_tmp.Rdata")

  # 3. Load file and extract information -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  load(paste(savedObj))

  if(!exists('DataBIOMOD') | !exists('Biomod.material')){
    stop("DataBIOMOD or Biomod.material Object not in your input object, please check it!")
  }

  ### tips to remove some compiling warnings
  if(!exists('Biomod.material')){ Biomod.material <- NULL }
  if(!exists('DataBIOMOD')){ DataBIOMOD <- NULL }
  if(!exists('Biomod.PA.sample')){ Biomod.PA.sample <- NULL }
  if(!exists('Evaluation.results.Roc')){ Evaluation.results.Roc <- NULL }
  if(!exists('Evaluation.results.TSS')){ Evaluation.results.TSS <- NULL }
  if(!exists('Evaluation.results.Kappa')){ Evaluation.results.Kappa <- NULL }
  if(!exists('VarImportance')){ VarImportance <- NULL }

  sp.names <- Biomod.material$species.names

  NewModelObj <- lapply(sp.names,function(sp.name){
    cat("\n\n", sp.name, 'run conversion...')

    dir.create(sp.name,showWarnings=FALSE)
    dir.create(paste(sp.name,"/.BIOMOD_DATA",sep=""),showWarnings=FALSE)

    models.out <- new('BIOMOD.models.out',
                      sp.name = sp.name,
                      expl.var.names = Biomod.material$VarNames,
                      rescal.all.models = FALSE)

    #   3.1 BIOMOD.formated.data creation
    cat("\n\tBIOMOD.formated.data creation")
    data <- BIOMOD.formated.data(sp = DataBIOMOD[,sp.name],
                                   env = DataBIOMOD[,Biomod.material$VarNames],
                                   xy = data.frame(),
                                   sp.name = sp.name)

    if(Biomod.material$NbRepPA > 0){
      PAtmp <- matrix(FALSE,
                      nrow=nrow(DataBIOMOD),
                      ncol=length(Biomod.PA.sample[[sp.name]]),
                      dimnames=list(NULL, names(Biomod.PA.sample[[sp.name]])))
      for(i in 1: ncol(PAtmp)){
        PAtmp[Biomod.PA.sample[[sp.name]][[i]],i] <- TRUE
      }
      data <- new('BIOMOD.formated.data.PA',
                    sp.name = data@sp.name,
                    coord = data@coord,
                    data.env.var = data@data.env.var,
                    data.species = data@data.species,
                    PA = data.frame(PAtmp))
    }
    # save Input Data
    save(data, file = paste(models.out@sp.name,"/.BIOMOD_DATA/formated.input.data",sep=""), compress=compress.arg)
    models.out@formated.input.data@inMemory <- FALSE
    models.out@formated.input.data@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/formated.input.data",sep="")
    rm(PAtmp,data)

    #   3.2 BIOMOD.Model.Options creation
    cat("\n\tBIOMOD.Model.Options creation")
    models.options <- BIOMOD_ModelingOptions()
    # save Model Options
    save(models.options, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.options",sep=""),  compress=compress.arg)
    models.out@models.options@inMemory <- FALSE
    models.out@models.options@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.options",sep="")
    rm(models.options)

    #   3.3 Converting Models Computed
    cat("\n\tConverting Models Computed")

    if(!file.exists(paste(path,"models/",sep=""))){
      stop("models directory doesn't exist!")
    }
    dir.create(paste(models.out@sp.name,"/models",sep=""), showWarnings=FALSE)
    if(file.exists(paste(path,"models/scaling_models",sep=""))){
      dir.create(paste(models.out@sp.name,"/models/scaling_models",sep=""), showWarnings=FALSE)
    }

    old.mod.computed <- list.files(path = paste(path,"models/",sep=""),
                               pattern = models.out@sp.name,
                               recursive = TRUE)

    new.mod.computed <- old.mod.computed
    new.mod.computed <- sapply(new.mod.computed, function(x){
      if(length(grep('PA',x)) > 0){ # pseudo absences done case
        if(length(grep('rep',x)) > 0){ # repetition
          x <- gsub('rep','RUN',x)
        } else{
          x <- paste(x,'_Full',sep='')
        }
        if(length(grep('Rmod_',x)) > 0){ # scaled models
          x <- paste(gsub('Rmod_','',x),'_scaled',sep='')

          x.str <- unlist(strsplit(gsub('scaling_models/','',x),'_'))
          x <- paste(x.str[1], x.str[3], x.str[4], x.str[2], x.str[length(x.str)], sep='_')
          x <- paste('scaling_models/',x,sep='')
        } else{
          x.str <- unlist(strsplit(x,'_'))
          x <- paste(x.str[1], x.str[3], x.str[4], x.str[2], sep='_')
        }
      } else{
        x <- gsub('_full','',x)
        if(length(grep('rep',x)) > 0){ # repetition
          x <- gsub('rep','RUN',x)
        } else{
          x <- paste(x,'_Full',sep='')
        }
        if(length(grep('Rmod_',x)) > 0){ # scaled models
          x <- paste(gsub('Rmod_','',x),'_scaled',sep='')

          x.str <- unlist(strsplit(gsub('scaling_models/','',x),'_'))
          x <- paste(x.str[1], x.str[3], x.str[2], x.str[length(x.str)], sep='_')
          x <- paste('scaling_models/',x,sep='')
        } else{
          x.str <- unlist(strsplit(x,'_'))
          x <- paste(x.str[1], x.str[3], x.str[2], sep='_')
        }
      }
      return(x)
    })


    # coping the files with appropriated names
    lapply(1:length(old.mod.computed), function(x){
      file.copy(from = paste(path, "models/", old.mod.computed[x], sep=""),
                to = paste(models.out@sp.name, "/models/", new.mod.computed[x], sep=""),
                overwrite = TRUE,
                recursive = FALSE,
                copy.mode = TRUE )
    })

    models.out@models.computed <- unique(as.character(gsub('_scaled','',
                                                    gsub('scaling_models/','',new.mod.computed))))
#     models.out@models.failed <- Biomod.material$calibration.failures

    #   3.3 Models evaluation conversion
    algo.choosen.id <- which( Biomod.material$algo.choice == TRUE)
    algo.choosen.names <- Biomod.material$algo[algo.choosen.id]

    if(Biomod.material$NbRunEval>0){
      run.eval.names <- c(paste('RUN',1:Biomod.material$NbRunEval, sep=''),'Full')
    } else{
      run.eval.names <- 'Full'
    }

    if(Biomod.material$NbRepPA>0){
      pa.data.names <- paste('PA',1:(Biomod.material$NbRun[which( Biomod.material$species.names == sp.name)] / (Biomod.material$NbRunEval+1) ),sep='')
    } else{
      pa.data.names <- 'AllData'
    }

    mod.eval.met <- c()
    if(!is.null(Evaluation.results.Roc)) mod.eval.met <- c(mod.eval.met, 'ROC')
    if(!is.null(Evaluation.results.TSS)) mod.eval.met <- c(mod.eval.met, 'TSS')
    if(!is.null(Evaluation.results.Kappa)) mod.eval.met <- c(mod.eval.met, 'KAPPA')

    models.evaluation <- array(data=NA,
                               dim=c(length(mod.eval.met), 4, length(algo.choosen.names),
                                     length(run.eval.names), length(pa.data.names)),
                               dimnames=list(mod.eval.met,
                                             c("Testing.data", "Cutoff", "Sensitivity", "Specificity"),
                                             algo.choosen.names,
                                             run.eval.names,
                                             pa.data.names))

#     outTmp <- lapply(pa.data.names, function(pdn){
    for( pdn in pa.data.names){
      pdnOld <- ifelse(pdn!='AllData',pdn,'full')
#       lapply(run.eval.names, function(ren){
      for(ren in run.eval.names){
        renOld <- ifelse(length(grep('RUN',ren)>0), gsub('RUN','rep',ren), '')
        runOldName <- ifelse(renOld!='',
                             paste(models.out@sp.name, pdnOld, renOld, sep='_'),
                             paste(models.out@sp.name, pdnOld, sep='_'))
        if(!is.null(Evaluation.results.Roc)){
          models.evaluation['ROC',,,ren,pdn] <- matrix(round(as.numeric(as.matrix(Evaluation.results.Roc[[runOldName]][,c('Cross.validation', 'Cutoff', "Sensitivity", "Specificity")])),digits=3), nrow=4, byrow=T)
        }

        if(!is.null(Evaluation.results.TSS)){
          models.evaluation['TSS',,,ren,pdn] <- matrix(round(as.numeric(as.matrix(Evaluation.results.TSS[[runOldName]][,c('Cross.validation', 'Cutoff', "Sensitivity", "Specificity")])),digits=3), nrow=4, byrow=T)
        }

        if(!is.null(Evaluation.results.Kappa)){
          models.evaluation['KAPPA',,,ren,pdn] <- matrix(round(as.numeric(as.matrix(Evaluation.results.Kappa[[runOldName]][,c('Cross.validation', 'Cutoff', "Sensitivity", "Specificity")]),digits=3)), nrow=4, byrow=T)
        }

      }#)
    }#)

    # save model evaluation
    save(models.evaluation, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.evaluation",sep=""),  compress=compress.arg)
    models.out@models.evaluation@inMemory <- TRUE
    models.out@models.evaluation@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.evaluation",sep="")
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)


    #   3.4 Models variable Importances

    # save model variables importances
    if(!is.null(VarImportance)){
      variables.importances <- t(VarImportance[[sp.name]])
      save(variables.importances, file = paste(models.out@sp.name,"/.BIOMOD_DATA/variables.importances",sep=""),  compress=compress.arg)
      models.out@variables.importances@inMemory <- TRUE
      models.out@variables.importances@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/variables.importances",sep="")
      models.out@variables.importances@val <- variables.importances
      rm(variables.importances)
    }


    #   3.5 Models Predictions
    if(!exists(paste('Pred_', sp.name,sep=''))){
      load(paste(path,'pred/Pred_',sp.name,sep=""))
    }
    if(exists(paste('Pred_', sp.name,sep=''))){

      eval(parse(text=paste('models.prediction <- Pred_', sp.name,sep='')))
#       eval(parse(text=paste('rm(Pred_', sp.name,')',sep='')))

      kept.dim <- dim(models.prediction)
      kept.dim[2] <- length(algo.choosen.names)
      kept.dimnames <- dimnames(models.prediction)
#       kept.dimnames[[1]] <- NULL
      kept.dimnames[[2]] <- algo.choosen.names

      models.prediction <- models.prediction[,algo.choosen.names,,]
      dim(models.prediction) <- kept.dim
      dimnames(models.prediction) <- kept.dimnames

      if(length(grep('PA',dimnames(models.prediction)[[4]] )) == 0){
        dimnames(models.prediction)[[4]] <- 'AllData'
      }
      if(dim(models.prediction)[3] > 1){
        vecTmp <- models.prediction[,,1,]
        models.prediction[,,1,] <- models.prediction[,,dim(models.prediction)[3],]
        models.prediction[,,dim(models.prediction)[3],] <- vecTmp
        rm(vecTmp)
        dimnames(models.prediction)[[3]] <- c(paste('RUN',1:(dim(models.prediction)[3] - 1),sep=''),'Full')
      } else{
        dimnames(models.prediction)[[3]] <- 'Full'
      }

      # save model predictions
      save(models.prediction, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction",sep=""),  compress=compress.arg)
      models.out@models.prediction@inMemory <- FALSE
      models.out@models.prediction@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction",sep="")
      rm(models.prediction)
    }

    return(models.out)

  })




  # xx. Reconstruct original workspace -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   rm(list=ls())
#   load("WS_tmp.Rdata")
#   file.remove("WS_tmp.Rdata")

  names(NewModelObj) <- Biomod.material$species.names

  cat("\n\n-=-=-=- Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
  return(NewModelObj)
}
##' @name BIOMOD_cv
##' @aliases BIOMOD_cv
##' 
##' @title Custom models cross-validation procedure
##' 
##' @description This function creates a DataSplitTable which could be used to evaluate models in Biomod with repeated
##'   k-fold cross-validation (cv) or stratified cv instead of repeated split sample runs
##' 
##' @param data            BIOMOD.formated.data object returned by BIOMOD_FormatingData
##' @param k               number of bins/partitions for k-fold cv
##' @param stratified.cv   logical. run a stratified cv 
##' @param stratify        stratification method of the cv. Could be "x", "y", "both" (default), "block" or the name of a predictor for environmental stratified cv.
##' @param balance         make balanced particions for "presences" (default) or "absences" (resp. pseudo-absences or background).
##' @param repetition      number of repetitions of k-fold cv (1 if stratified.cv=TRUE)
##' @param do.full.models  if true, models calibrated and evaluated with the whole dataset are done
##' 
##' @details
##'   Stratified cv could be used to test for model overfitting and for assessing transferability in geographic and environmental space. 
##'   If balance = "presences" presences are divided (balanced) equally over the particions (e.g. Fig. 1b in Muscarelly et al. 2014).
##'   Pseudo-Absences will however be unbalanced over the particions especially if the presences are clumped on an edge of the study area.
##'   If balance = "absences" absences (resp. Pseudo-Absences or background) are divided (balanced) as equally as possible for the particions
##'   (geographical balanced bins given that absences are spread over the study area equally, approach similar to Fig. 1 in Wenger et Olden 2012).
##'   Presences will however be unbalanced over the particians. Be careful: If the presences are clumped on an edge of the study area it is possible that all presences are in one bin.
##' 
##' @return
##' DataSplitTable matrix with k*repetition (+ 1 for Full models if  do.full.models = TRUE) columns for BIOMOD_Modeling function.
##' Stratification "x" and "y" was described in Wenger and Olden 2012. While Stratification "y" uses k partitions along the y-gradient, "x" does the same for the x-gradient and "both" combines them.
##' Stratification "block" was described in Muscarella et al. 2014. For bins of equal number are partitioned (bottom-left, bottom-right, top-left and top-right).
##' 
##' @references
##' Muscarella, R., Galante, P.J., Soley-Guardia, M., Boria, R.A., Kass, J.M., Uriarte, M. & Anderson, R.P. (2014). ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, \bold{5}, 1198-1205.
##' Wenger, S.J. & Olden, J.D. (2012). Assessing transferability of ecological models: an underappreciated aspect of statistical validation. \emph{Methods in Ecology and Evolution}, \bold{3}, 260-267.
##' 
##' @author Frank Breiner \email{frank.breiner@wsl.ch}
##' 
##' @seealso
##' \code{\link[ENMeval]{get.block}}
##'
##' @examples
##' \dontrun{
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"))
##' head(DataSpecies)
##' 
##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio7.grd", 
##'                              package="biomod2"),  
##'                 system.file( "external/bioclim/current/bio11.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio12.grd", 
##'                              package="biomod2"))
##' 
##' # 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' 
##' # 3. Creating DataSplitTable
##' 
##' DataSplitTable <- BIOMOD_cv(myBiomodData, k=5, rep=2, do.full.models=F)
##' DataSplitTable.y <- BIOMOD_cv(myBiomodData,stratified.cv=T, stratify="y", k=2)
##' colnames(DataSplitTable.y)[1:2] <- c("RUN11","RUN12")
##' DataSplitTable <- cbind(DataSplitTable,DataSplitTable.y)
##' head(DataSplitTable)
##' 
##' # 4. Doing Modelisation
##' 
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##'                                      models = c('RF'), 
##'                                      models.options = myBiomodOption, 
##'                                      DataSplitTable = DataSplitTable,
##'                                      VarImport=0, 
##'                                      models.eval.meth = c('ROC'),
##'                                      do.full.models=FALSE,
##'                                      modeling.id="test")
##' 
##' ## get cv evaluations
##' eval <- get_evaluations(myBiomodModelOut,as.data.frame=T)
##' 
##' eval$strat <- NA
##' eval$strat[grepl("13",eval$Model.name)] <- "Full"
##' eval$strat[!(grepl("11",eval$Model.name)|
##'              grepl("12",eval$Model.name)|
##'              grepl("13",eval$Model.name))] <- "Random"
##' eval$strat[grepl("11",eval$Model.name)|grepl("12",eval$Model.name)] <- "Strat"
##' 
##' boxplot(eval$Testing.data~ eval$strat, ylab="ROC AUC")
##' }


BIOMOD_cv <-
  function (data, k = 5, repetition = 5, do.full.models = TRUE, 
            stratified.cv = FALSE, stratify = "both", balance = "pres") 
  {
    DataSplitTable.y <- DataSplitTable.x <- DataSplitTable <- NULL
    if (stratified.cv) {
      repetition <- 1
      if (balance == "absences") {
        balance <- data@data.species == 1 | data@data.species == 
          0
      }
      else {
        balance <- data@data.species == 1
      }
      if (stratify == "x" | stratify == "both") {
        DataSplitTable.x <- matrix(NA, nrow(data@coord), 
                                   k)
        bands <- quantile(data@coord[balance, 1], probs = seq(0, 
                                                              100, 100/k)/100)
        bands[1] <- bands[1] - 1
        bands[k + 1] <- bands[k + 1] + 1
        for (i in 1:k) {
          DataSplitTable.x[, i] <- data@coord[, 1] >= bands[i] & 
            data@coord[, 1] < bands[i + 1]
        }
        if (stratify == "x") {
          DataSplitTable <- DataSplitTable.x
        }
      }
      if (stratify == "y" | stratify == "both") {
        DataSplitTable.y <- matrix(NA, nrow(data@coord), 
                                   k)
        bands <- quantile(data@coord[balance, 2], probs = seq(0, 
                                                              100, 100/k)/100)
        bands[1] <- bands[1] - 1
        bands[k + 1] <- bands[k + 1] + 1
        for (i in 1:k) {
          DataSplitTable.y[, i] <- data@coord[, 2] >= bands[i] & 
            data@coord[, 2] < bands[i + 1]
        }
        if (stratify == "y") {
          DataSplitTable <- DataSplitTable.y
        }
      }
      if (stratify == "both") {
        DataSplitTable <- cbind(DataSplitTable.x, DataSplitTable.y)
      }
    if (stratify == "block") {
      DataSplitTable <- as.data.frame(matrix(NA, nrow(data@coord), 
                                             4))

      blocks<-ENMeval::get.block(data@coord[data@data.species==1,],
                                 data@coord[data@data.species==0,])
      
      for(i in 1:4){
        DataSplitTable[data@data.species == 1,i] <-  blocks[[1]]!=i     
        DataSplitTable[data@data.species == 0,i] <-  blocks[[2]]!=i     
      }
    }      
      if (stratify != "block" & stratify != "x" & stratify != 
          "y" & stratify != "both") {
        DataSplitTable2 <- as.data.frame(matrix(NA, nrow(data@coord), k))
        bands <- quantile(data@data.env.var[balance, stratify], 
                          probs = seq(0, 100, 100/k)/100)
        bands[1] <- bands[1] - 1
        bands[k + 1] <- bands[k + 1] + 1
        for (i in 1:k) {
          DataSplitTable2[, i] <- data@data.env.var[balance, 
                                                    stratify] <= bands[i] | data@data.env.var[balance, 
                                                                                              stratify] > bands[i + 1]
        }
      }
    }
    else {
      for (rep in 1:repetition) {
        fold <- dismo::kfold(data@data.species, by = data@data.species, 
                             k = k)
        for (i in 1:k) {
          DataSplitTable <- cbind(DataSplitTable, fold != 
                                    i)
        }
      }
    }
    if(stratify != "block"){
      colnames(DataSplitTable) <- paste("RUN", 1:(k * repetition), 
                                        sep = "")
      if (do.full.models == TRUE) {
        DataSplitTable <- cbind(DataSplitTable, T)
        colnames(DataSplitTable)[k * repetition + 1] <- "Full"
      }
    }else{
      colnames(DataSplitTable) <- paste("RUN", 1:4, 
                                        sep = "")    
      if (do.full.models == TRUE) {
        DataSplitTable <- cbind(DataSplitTable, T)
        colnames(DataSplitTable)[5] <- "Full"
      }
    }
    
    return(DataSplitTable)
  }
'BIOMOD_EnsembleForecasting' <- function( EM.output,
                                          projection.output = NULL,
                                          new.env = NULL,
                                          xy.new.env = NULL,
                                          selected.models = 'all',
                                          proj.name = NULL,
                                          binary.meth = NULL,
                                          filtered.meth = NULL,
                                          compress = TRUE,
                                          ...){
  .bmCat("Do Ensemble Models Projections")
  
  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  args <- list(...)
  
  args_checked <- .BIOMOD_EnsembleForecasting.check.args( EM.output,
                                                          projection.output,
                                                          new.env,
                                                          selected.models,
                                                          proj.name,
                                                          total.consensus=FALSE,
                                                          binary.meth,
                                                          filtered.meth )
  
#   total.consensus <- args_checked$total.consensus
  proj.name <- args_checked$proj.name
  selected.models <- args_checked$selected.models
    
  output.format <- args$output.format # raster output format
  compress <- args$compress # compress or not output
  do.stack <- args$do.stack # save raster as stack or layers
  keep.in.memory <- args$keep.in.memory # store results on memory or only on hard drive
  on_0_1000 <- args$on_0_1000 # convert 0-1 predictions on 0-1000 scale to limit memory consuming


  
  
  if(is.null(output.format)){
    if(length(projection.output)){
      if(projection.output@type != 'RasterStack')
        output.format <- ".RData"
      else
        output.format <- ".grd"
    } else{
      if(! inherits(new.env,'Raster'))
        output.format <- ".RData"
      else
        output.format <- ".grd"
    }

  }
  
  if(is.null(compress)) compress <- FALSE
  
  if(is.null(do.stack)) {
    do.stack <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than projection.output ones
    if(!is.null(projection.output)){
      if(all(grepl("individual_projections", projection.output@proj@link))){
        do.stack <- FALSE
      }
    }
  }
  
  if(is.null(keep.in.memory)){
    keep.in.memory <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than projection.output ones
    if(!is.null(projection.output)){
      keep.in.memory <- projection.output@proj@inMemory
    }
  } 
  
  if(is.null(xy.new.env)) { 
    if(!is.null(projection.output)){
      xy.new.env <- projection.output@xy.coord
    } else{
      xy.new.env <- matrix()
    }
  } 

  if(is.null(on_0_1000)) on_0_1000 = TRUE # by default outputs are return on a 0 - 1000 scale 
  
  rm(list=c('args_checked','args'))
  
  # get argument from projection
  proj_out <- new('BIOMOD.projection.out',
                  proj.names = proj.name,
                  sp.name =  EM.output@sp.name,
                  expl.var.names = EM.output@expl.var.names,
                  models.projected = selected.models,
#                   scaled.models = modeling.output@rescal.all.models,
                  xy.coord = xy.new.env,
                  modeling.object.id = EM.output@modeling.id)
  
  proj_out@modeling.object@link = EM.output@link
  
  proj_is_raster <- FALSE
  if(inherits(new.env, 'Raster')){
    proj_is_raster <- TRUE
  } else if( length(projection.output) ){
    if(inherits(projection.output@proj, 'BIOMOD.stored.raster.stack')){
      proj_is_raster <- TRUE
    }
  }
  
  
  if(proj_is_raster){
    proj_out@proj <- new('BIOMOD.stored.raster.stack')
  } else{
    proj_out@proj <- new('BIOMOD.stored.array')
    do.stack = TRUE
  }
  
  # 1.c creating output directory
  dir.create(file.path(EM.output@sp.name,paste("proj_", proj.name, sep="")), 
             showWarnings = FALSE, recursive = TRUE, mode = "777")
  
  indiv_proj_dir <- file.path(EM.output@sp.name,paste("proj_", proj.name, sep=""), "individual_projections")
#   if(!do.stack){
#     rasterOptions(todisk=T)
    dir.create(indiv_proj_dir, 
               showWarnings = FALSE, recursive = TRUE, mode = "777")
#   }
  
  saved.files <- c()
                                                  
  
  ### get needed models prediction ###
  needed_predictions <- get_needed_models(EM.output, selected.models=selected.models)

  if (length(projection.output)){
    formal_pred <- get_predictions(projection.output, full.name=needed_predictions, as.data.frame=ifelse(projection.output@type=='array',T,F) )
  } else{
    # make prediction according to given environment
    tmp_dir <- paste('Tmp', format(Sys.time(), "%s"), sep="")
    formal_pred <- BIOMOD_Projection( modeling.output = load_stored_object(EM.output@models.out.obj),
                                      new.env = new.env,
                                      proj.name = tmp_dir,
                                      xy.new.env = NULL,
                                      selected.models = needed_predictions,
                                      compress = TRUE,
                                      build.clamping.mask = F,
                                      do.stack=T, silent = T, on_0_1000 = on_0_1000 )
    # getting the results
    formal_pred <- get_predictions(formal_pred, full.name=needed_predictions, as.data.frame=ifelse(inherits(new.env,'Raster'),F,T))
    
    # remove tmp directory
    unlink(file.path(EM.output@sp.name,paste("proj_",tmp_dir,sep="")),recursive = TRUE, force = TRUE)
  }
  
  ef.out <- NULL
  # 2. Do the ensemble modeling
  for( em.comp in EM.output@em.computed[which(EM.output@em.computed %in% selected.models)]){
    cat("\n\n\t> Projecting", em.comp, "...")
    model.tmp <- NULL
    file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp,output.format))
    BIOMOD_LoadModels(EM.output, full.name=em.comp, as='model.tmp')
    if(inherits(formal_pred,'Raster')){
      ef.tmp <- predict(model.tmp, 
                        formal_predictions = raster::subset(formal_pred, subset=model.tmp@model, drop=FALSE),
                        on_0_1000 = on_0_1000, 
                        filename = ifelse(output.format == '.RData', '', file_name_tmp))
    } else {
      ef.tmp <- predict(model.tmp, formal_predictions = formal_pred[,model.tmp@model, drop=FALSE], on_0_1000 = on_0_1000)
    }
        
    if(inherits(ef.tmp,'Raster')){
      if(do.stack){
        if(length(ef.out)) ef.out <- stack(ef.out,ef.tmp) else ef.out <- raster::stack(ef.tmp)
      } else {
        file_name_tmp <- file.path(indiv_proj_dir,paste(em.comp,output.format,sep=""))
        if(output.format== '.RData'){
          save(ef.tmp, file=file_name_tmp, compress=compress)
        } 
        saved.files <- c(saved.files, file_name_tmp)
      }
    } else{
      ef.out <- cbind(ef.out,ef.tmp)
    } 
  }
  
  proj_out@models.projected <- EM.output@em.computed[which(EM.output@em.computed %in% selected.models)]
  
  if(do.stack){
    if( inherits(ef.out, "Raster") ) {
      names(ef.out) <- proj_out@models.projected
    } else {
      colnames(ef.out) <- proj_out@models.projected
    }
    # save object
    file_name_tmp <- file.path(EM.output@sp.name,paste("proj_", proj.name, sep=""),paste("proj_", proj.name,"_",EM.output@sp.name,"_ensemble",output.format,sep=""))
    if(output.format== '.RData'){
      save(ef.out, file=file_name_tmp, compress=compress)
    } else if( inherits(ef.out, "Raster") ){
      ## TODO : define the raster dataformat (depends if em.cv has been computed)
      writeRaster(ef.out,filename=file_name_tmp, overwrite=TRUE)
    }
    saved.files <- c(saved.files, file_name_tmp)
    proj_out@proj@link <- file_name_tmp
  } else {
    proj_out@proj@link <- saved.files #EM.output@em.computed
  }
  
  if(!is.null(ef.out)){
    proj_out@proj@val <- ef.out
    proj_out@proj@inMemory <- TRUE
  }
  
  
  
  ## binary && filtering transformations #############################################################
  if(length(binary.meth) | length(filtered.meth)){
    cat("\n")
    eval.meth <- unique(c(binary.meth,filtered.meth))
    
    ## get all treshold
    thresholds <- sapply(selected.models, function(x){
      get_evaluations(EM.output)[[x]][eval.meth, "Cutoff"]  
    })
    
    ## convert thresholds deopending on the chosen scale
    if(! on_0_1000 ) { thresholds <- thresholds / 1000 }
    
    ## do binary transformation
    for(eval.meth in binary.meth){
      cat("\n\t> Building", eval.meth,"binaries")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          proj_bin <- BinaryTransformation(raster(file.tmp, RAT=FALSE),thres.tmp)
          writeRaster(x = proj_bin,
                      filename = sub(output.format, paste("_",eval.meth,"bin", output.format, sep=""), file.tmp), 
                      overwrite=TRUE ) ## , datatype = "INT2S",NAflag=-9999)
        }
      } else {
        assign(x = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep=""),
               value = BinaryTransformation(ef.out,thresholds))
        
        if(output.format == '.RData'){
          save(list = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep=""), 
               file = file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", output.format ,sep="")), compress=compress)   
        } else {
          writeRaster(x=get(paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep="")),
                      filename=file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", output.format ,sep="")), 
                      overwrite=TRUE) ## , datatype = "INT2S", NAflag=-9999)
                      
        }
        
        rm(list=paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep=""))
      }
    }
    
    ## do filtered transformation
    for(eval.meth in filtered.meth){
      cat("\n\t> Building", eval.meth,"filtered")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          ## TODO : define the raster dataformat (depends if em.cv has been computed)
          filt_proj <- FilteringTransformation(raster(file.tmp, RAT=FALSE),thres.tmp)
          writeRaster(x = filt_proj,
                      filename = sub(output.format, paste("_",eval.meth,"filt", output.format, sep=""), file.tmp), 
                      overwrite=TRUE)
        }
      } else {
        assign(x = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep=""),
               value = FilteringTransformation(ef.out,thresholds))
        
        if(output.format == '.RData'){
          save(list = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep=""), 
               file = file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", output.format ,sep="")), compress=compress)   
        } else {
          ## TODO : define the raster dataformat (depends if em.cv has been computed)
          writeRaster(x=get(paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep="")),
                      filename=file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", output.format ,sep="")), overwrite=TRUE)
        }
        
        rm(list=paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep=""))
      }
    }
  }
  
  
  # save object
  if(!keep.in.memory){
    proj_out <- free(proj_out)
  }
  
  assign(paste(EM.output@sp.name,".", proj.name, ".ensemble.projection.out", sep=""), proj_out)
  save(list = paste(EM.output@sp.name,".", proj.name, ".ensemble.projection.out", sep=""),
       file = file.path(EM.output@sp.name, paste("proj_", proj.name, sep=""), paste(EM.output@sp.name,".", proj.name, ".ensemble.projection.out", sep="")))
  
  cat("\n")
  .bmCat("Done")
  return(proj_out)
}







# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.BIOMOD_EnsembleForecasting.check.args <- function( EM.output,
                                                    projection.output,
                                                    new.env,
                                                    selected.models,
                                                    proj.name,
                                                    total.consensus,
                                                    binary.meth,
                                                    filtered.meth ){

  # check needed data are referred
  if( (is.null(projection.output) & is.null(new.env)) | (!is.null(projection.output) & !is.null(new.env))){
    stop("You have to refer at one of 'projection.output' or 'new.env' argument")
  }
  
  if(!inherits(EM.output,"BIOMOD.EnsembleModeling.out")){
    stop("EM.output must be a 'BIOMOD.EnsembleModeling.out' object")
  }
  
  # check all needed predictions are available
  needed_pred <- get_needed_models(EM.output, selected.models=selected.models)  
  
  if(!is.null(projection.output)){
    if(!inherits(projection.output, "BIOMOD.projection.out")){
      stop("projection.output must be a 'BIOMOD.projection.out' object")
    }
    missing_pred <- needed_pred[! (needed_pred %in% projection.output@models.projected)]
    if( length(missing_pred) ){
      stop("Some models prediction missing :", toString(missing_pred))
    }
  }
  
  ## selected.models
  if(selected.models[1] == 'all'){
    selected.models <- get_built_models(EM.output)
  } else{
    selected.models <- intersect(selected.models, get_built_models(EM.output))
  }
  if(length(selected.models) < 1){
    stop('No models selected')
  }
  
  ## projection name
  if(!length(proj.name) & !length(projection.output)){
    stop("You have to give a valid 'proj.name' if you don't work with projection.output")
  } else if(!length(proj.name)){
    proj.name <- projection.output@proj.names
  }
  
  
                              
  if(total.consensus){
    if(length(EM.output@em.computed) < 2){
      cat("\n      ! Total consensus projection was switched off because only one Ensemble modeling was done")
      total.consensus <- FALSE
    }
  }
    
  if(!is.null(binary.meth)){
#     if(sum(!(binary.meth %in% EM.output@eval.metric))){
#       stop(paste("binary methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
#                  toString(EM.output@eval.metric)," )", sep=""))
#     }
  }
  
  if(!is.null(filtered.meth)){
#     if(sum(!(filtered.meth %in% EM.output@eval.metric))){
#       stop(paste("filtering methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
#                  toString(EM.output@eval.metric)," )", sep=""))
#     }
  }
                                                    
  return(list(projection.output = projection.output,
              EM.output = EM.output,
              selected.models = selected.models,
              proj.name = proj.name,
              total.consensus = total.consensus,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth))
  
}
'BIOMOD_EnsembleModeling' <- function( modeling.output,
                                       chosen.models = 'all',
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = 'all',
                                       eval.metric.quality.threshold = NULL,
                                       models.eval.meth = c('KAPPA','TSS','ROC'),
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = FALSE,
                                       prob.mean.weight.decay = 'proportional',
                                       VarImport=0){
  .bmCat("Build Ensemble Models")
  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  args <- .BIOMOD_EnsembleModeling.check.args( modeling.output,
                                               chosen.models,
                                               eval.metric,
                                               eval.metric.quality.threshold,
                                               models.eval.meth,
                                               prob.mean,
                                               prob.cv,
                                               prob.ci,
                                               prob.ci.alpha,
                                               prob.median,
                                               committee.averaging,
                                               prob.mean.weight,
                                               prob.mean.weight.decay,
                                               em.by)

  modeling.output <- args$modeling.output
  chosen.models <- args$chosen.models
  eval.metric <- args$eval.metric
  eval.metric.quality.threshold <- args$eval.metric.quality.threshold
  models.eval.meth <- args$models.eval.meth
  prob.mean <- args$prob.mean
  prob.cv <- args$prob.cv
  prob.ci <- args$prob.ci
  prob.ci.alpha <- args$prob.ci.alpha
  prob.median <- args$prob.median
  committee.averaging <- args$committee.averaging
  prob.mean.weight <- args$prob.mean.weight
  prob.mean.weight.decay  <- args$prob.mean.weight.decay
  em.by <- args$em.by

  rm('args')

  em.avail <- c('prob.mean', 'prob.cv', 'prob.ci.inf', 'prob.ci.sup', 'prob.median', 'committee.averaging', 'prob.mean.weight')
  em.algo <- em.avail[c(prob.mean, prob.cv, prob.ci, prob.ci, prob.median, committee.averaging, prob.mean.weight)]

  # create a EM option list
  Options <- list(em.by=em.by)
  expl_var_type = get_var_type(get_formal_data(modeling.output,'expl.var'))
  expl_var_range = get_var_range(get_formal_data(modeling.output,'expl.var'))


  # 1b. creating output object and begin to fill it
#   EM <- list()
  EM <- new('BIOMOD.EnsembleModeling.out',
            sp.name = modeling.output@sp.name,
            expl.var.names = modeling.output@expl.var.names,
            em.by = em.by,
            modeling.id = modeling.output@modeling.id
#             models.out.obj = new('BIOMOD.stored.models.out',
#                                  inMemory = FALSE,
#                                  link = paste(modeling.output@sp.name,"/",modeling.output@sp.name,".models.out",sep="")),
#             eval.metric = eval.metric,
#             eval.metric.quality.threshold = eval.metric.quality.threshold#,
#             em.ci.alpha = prob.ci.alpha
            )

  EM@models.out.obj@link <- file.path(modeling.output@sp.name,paste(modeling.output@sp.name,".", modeling.output@modeling.id,".models.out",sep="") )

  # 2. doing Ensemble modeling

  ## 2.1 make a list of models names that will be combined together according to by argument.
  em.mod.assemb <- .em.models.assembling(chosen.models, em.by)

  for(assemb in names(em.mod.assemb) ){
    cat("\n\n  >", assemb, "ensemble modeling")
    models.kept <- em.mod.assemb[[assemb]]

    #### defined data that will be used for models performances calculation ####
    if(modeling.output@has.evaluation.data){
      eval.obs <- get_formal_data(modeling.output,'eval.resp.var')
      eval.expl <- get_formal_data(modeling.output,'eval.expl.var')
    }

    ##### !!!!!! TO DO -> select appropriate part of dataset according to em.by
    obs <-  get_formal_data(modeling.output,'resp.var')
    expl <- get_formal_data(modeling.output,'expl.var')

    ## subselection of observations according to dataset used to produce ensemble models
    if(em.by %in% c("PA_dataset",'PA_dataset+algo','PA_dataset+repet')){
      if(unlist(strsplit(assemb,"_"))[3] != 'AllData'){
        if(inherits(get_formal_data(modeling.output), "BIOMOD.formated.data.PA")){
          kept_cells <- get_formal_data(modeling.output)@PA[, unlist(strsplit(assemb,"_"))[3]]
        } else {
          kept_cells <- rep(T, length(obs))
        }

        obs <- obs[kept_cells]
        expl <- expl[kept_cells, ,drop=F]
      }
    }

    ## subselection of observations according to dataset used to produce ensemble models is done at evaluation step


#     if(em.by %in% c("algo","all") ){
#       ## we need to take all data even if it is much better to have
#       obs <-  get_formal_data(modeling.output,'resp.var')
#       expl <- get_formal_data(modeling.output,'expl.var')
#     }
    # remove na
    obs[is.na(obs)] <- 0




    #### get needed models predictions ####
    needed_predictions <- .get_needed_predictions(modeling.output, em.by, models.kept, eval.metric, eval.metric.quality.threshold)
    # if no prediction selected => swith to next model
    if(!length(needed_predictions)) next

    ## loop on evaluation metrics ##
    for(eval.m in eval.metric){

      # define model name
#       base_model_name <- paste(modeling.output@sp.name,"_",assemb,"_",eval.m,'_' ,sep="")
      models.kept <- needed_predictions$models.kept[[eval.m]]
      models.kept.scores <- needed_predictions$models.kept.scores[[eval.m]]

      ## Loop over em.algo ##

      for(algo in em.algo){
        #### Models building ####

        # 1. Mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.mean'){
          cat("\n   > Mean of probabilities...")

#           model_name <- paste(base_model_name,"EMmean",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMmeanBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMmean_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMmean',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id)
        }

        # 2. CV of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.cv'){
          cat("\n   > Coef of variation of probabilities...")
#           model_name <- paste(base_model_name,"EMcv",sep="")
          model_name <-paste(modeling.output@sp.name,"_","EMcvBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMcv_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMcv',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id)
        }

        # 3. Median of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.median'){
          cat("\n   > Median of probabilities...")
#           model_name <- paste(base_model_name,"EMmedian",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMmedianBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMmedian_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMmedian',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id)
        }

        # 4. CI of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.ci.inf'){
          cat("\n   > Confidence Interval...")

          ## Quantile inferior
#           model_name <- paste(base_model_name,"EMciInf",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMciInfBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMci_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMci',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           alpha = prob.ci.alpha,
                           side = 'inferior')
        }

        if(algo == 'prob.ci.sup'){
          ## Quantile superior
#           model_name <- paste(base_model_name,"EMciSup",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMciSupBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMci_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMci',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           alpha = prob.ci.alpha,
                           side = 'superior')

        }

        # 5. Comitee averaging of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(algo == 'committee.averaging'){
          cat("\n   >  Committee averaging...")
#           model_name <- paste(base_model_name,"EMca",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMcaBy",eval.m, "_", assemb ,sep="")

          models.kept.tresh <- unlist(lapply(models.kept, function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(get_evaluations(modeling.output)[eval.m, "Cutoff", mod, run, dat])
          }))
          names(models.kept.tresh) <- models.kept

          ## remove models if some thresholds are undefined
          to_keep <- is.finite(models.kept.tresh)

          model.bm <- new("EMca_biomod2_model",
                           model = models.kept[to_keep],
                           model_name = model_name,
                           model_class = 'EMca',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           tresholds = models.kept.tresh[to_keep])

        }

        # 6. weighted mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(algo == 'prob.mean.weight'){
          cat("\n   > Probabilities weighting mean...")
#           model_name <- paste(base_model_name,"EMwmean",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMwmeanBy",eval.m, "_", assemb ,sep="")

          # remove SRE models if ROC
          models.kept.tmp <- models.kept
          models.kept.scores.tmp <- models.kept.scores
#           prediction.kept.tmp <- prediction.kept

          if(eval.m == 'ROC'){
            sre.id <- grep("_SRE", models.kept)
            if(length(sre.id)>0){
              cat("\n      ! SRE modeling were switched off")
              models.kept.tmp <- models.kept[-sre.id]
              models.kept.scores.tmp <- models.kept.scores[-sre.id]
#               prediction.kept.tmp <- prediction.kept[,models.kept]
            }
          }

          ## remove models if score is not defined
          models.kept.tmp <- models.kept.tmp[is.finite(models.kept.scores.tmp)]
          models.kept.scores.tmp <- models.kept.scores.tmp[is.finite(models.kept.scores.tmp)]

          # weights are "decay" times decreased for each subsequent model in model quality order.
          models.kept.scores.tmp <- round(models.kept.scores.tmp, 3) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.

          # dealing with numerical decay
          cat("\n\t\t", " original models scores = ", models.kept.scores.tmp)
          if(is.numeric(prob.mean.weight.decay)){
            DecayCount <- sum(models.kept.scores.tmp>0)
            WOrder <- order(models.kept.scores.tmp, decreasing=T)
            Dweights <- models.kept.scores.tmp
            ## old version
            # for(J in 1:DecayCount) Dweights[WOrder[J]] <- (DecayCount - J + 1) * prob.mean.weight.decay
            ## end old version
            for(J in 1:DecayCount) Dweights[WOrder[J]] <- I(prob.mean.weight.decay^(DecayCount - J + 1))
            #If 2 or more score are identical -> make a mean weight between the ones concerned
            for(J in 1:length(models.kept.scores.tmp)){
              if(sum(models.kept.scores.tmp[J]==models.kept.scores.tmp)>1) Dweights[which(models.kept.scores.tmp[J]==models.kept.scores.tmp)] <- mean(Dweights[which(models.kept.scores.tmp[J]==models.kept.scores.tmp)])
            }
            models.kept.scores.tmp <- round(Dweights, digits=3)
            rm(list=c('Dweights','DecayCount','WOrder'))
          } else if ( is.function(prob.mean.weight.decay) ){ # dealing with function decay
            models.kept.scores.tmp <- sapply(models.kept.scores.tmp, prob.mean.weight.decay)
          }

          ### Standardise model weights
          models.kept.scores.tmp <- round(models.kept.scores.tmp/sum(models.kept.scores.tmp, na.rm=T), digits=3)

          cat("\n\t\t", " final models weights = ", models.kept.scores.tmp)

          model.bm <- new("EMwmean_biomod2_model",
                           model = models.kept.tmp,
                           model_name = model_name,
                           model_class = 'EMwmean',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           penalization_scores = models.kept.scores.tmp)
        }

        #### Models Evaluation ####
        pred.bm <- predict(model.bm, expl, formal_predictions=needed_predictions$predictions[,model.bm@model, drop=F], on_0_1000 = T )

        ## store models prediction on the hard drive ---------------------------
        ## create the suitable directory architecture
        pred.bm.name <- paste0(model_name, ".predictions")
        pred.bm.outfile <- file.path(model.bm@resp_name, ".BIOMOD_DATA", model.bm@modeling.id,
                                     "ensemble.models", "ensemble.models.predictions",
                                     pred.bm.name)
        dir.create(dirname(pred.bm.outfile), showWarnings = FALSE, recursive = TRUE)
        ## save models predictions
        assign(pred.bm.name, pred.bm)
        save(list = pred.bm.name, file = pred.bm.outfile, compress = TRUE)
        rm(list = pred.bm.name)
        ## end strore models preciciton on the hard drive ----------------------

        if(exists('eval.obs') & exists('eval.expl')){
          eval_pred.bm <- predict(model.bm, eval.expl)
          ## store models prediction on the hard drive -------------------------
          pred.bm.name <- paste0(model_name, ".predictionsEval")
          ## save models predictions
          assign(pred.bm.name, eval_pred.bm)
          save(list = pred.bm.name, file = pred.bm.outfile, compress = TRUE)
          rm(list = pred.bm.name)
          ## end strore models preciciton on the hard drive --------------------
        }


        # Model evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if( length(models.eval.meth) ){
          cat("\n\t\t\tEvaluating Model stuff...")

          if(algo == 'prob.cv'){ ## switch of evalutaion process
            cross.validation <- matrix(NA,4,length(models.eval.meth),
                                       dimnames = list(c("Testing.data","Cutoff","Sensitivity", "Specificity"),
                                                       models.eval.meth))
          } else {
            if(em.by == "PA_dataset+repet"){ ## select the same evaluation data than formal models
              ## get formal models calib/eval lines
              calib_lines <- get_calib_lines(modeling.output)
              ## get info on wich dataset and which repet this ensemble model is based on
              pa_dataset_id <- paste("_", unlist(strsplit(assemb,"_"))[3], sep="")
              repet_id <- paste("_", unlist(strsplit(assemb,"_"))[2], sep="")
              ## define and extract the subset of points model will be evaluated on
              if (repet_id == "_Full"){
                eval_lines <- rep(T, length(pred.bm))
              } else {
                eval_lines <- ! na.omit(calib_lines[ , repet_id, pa_dataset_id])
                ## trick to detect when it is a full model but with a non common name
                if(all(!eval_lines)){ ## i.e. all lines used for calib => full model
                  eval_lines <- !eval_lines
                }
              }
            } else {
              eval_lines <- rep(T, length(pred.bm))
            }

            cross.validation <- sapply(models.eval.meth,
                                       Find.Optim.Stat,
                                       Fit = pred.bm[eval_lines],
                                       Obs = obs[eval_lines])
            rownames(cross.validation) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
          }



          if(exists('eval_pred.bm')){

            if(algo == 'prob.cv'){ ## switch of evalutaion process
              cross.validation <- matrix(NA,5,length(models.eval.meth),
                                         dimnames = list(c("Testing.data","Evaluating.data","Cutoff","Sensitivity", "Specificity"),
                                                         models.eval.meth))
            } else {
              true.evaluation <- sapply(models.eval.meth, function(x){
                                        Find.Optim.Stat(
                                        Fit = eval_pred.bm * 1000,
                                        Obs = eval.obs,
                                        Fixed.thresh = cross.validation["Cutoff",x])})

              cross.validation <- rbind(cross.validation["Testing.data",], true.evaluation)
              rownames(cross.validation) <- c("Testing.data","Evaluating.data","Cutoff","Sensitivity", "Specificity")
            }
          }

          ## store results
          cross.validation <- t(round(cross.validation,digits=3))
          model.bm@model_evaluation <- cross.validation

          ## remove useless objects
          rm(list=c('cross.validation') )
        }
        # End evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

        #### Var Importance calculation ####
        if (VarImport > 0){ # do Varimp stuff
          cat("\n\t\t\tEvaluating Predictor Contributions...", "\n")
          variables.importance <- variables_importance(model.bm, expl, nb_rand=VarImport)
          model.bm@model_variables_importance <- variables.importance$mat
          ## remove useless objects
          rm(list=c('variables.importance') )
        }

        #### Models saving #####
        assign(model_name,model.bm)
        save(list=model_name,file=file.path(modeling.output@sp.name,
                                            "models",
                                            modeling.output@modeling.id,
                                            model_name))

        #### Add to sumary objects ####
        EM@em.models <- c(EM@em.models, model.bm)
        EM@em.computed <- c(EM@em.computed, model_name)


      }
    }
  }

  ### fix models names ###
  names(EM@em.models) <- EM@em.computed

  model.name <- paste(EM@sp.name, '.', EM@modeling.id, 'ensemble.models.out', sep="")
  assign(x=model.name,
         value=EM)
  save(list=model.name,
       file=file.path(EM@sp.name,model.name))

  .bmCat("Done")
  return(EM)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.BIOMOD_EnsembleModeling.check.args <- function(  modeling.output,
                                                   chosen.models,
                                                   eval.metric,
                                                   eval.metric.quality.threshold,
                                                   models.eval.meth,
                                                   prob.mean,
                                                   prob.cv,
                                                   prob.ci,
                                                   prob.ci.alpha,
                                                   prob.median,
                                                   committee.averaging,
                                                   prob.mean.weight,
                                                   prob.mean.weight.decay,
                                                   em.by ){
  # 1. modeling.output checking
  if(!(inherits(modeling.output, "BIOMOD.models.out"))){
    stop("Invalid modeling.output argument !\nIt must be a 'BIOMOD.models.out' object")
  }

  # 2. chosen.models checking
  if(!length(chosen.models) | (length(chosen.models)==1 & chosen.models[1] == 'all')){ # select all models
    cat("\n   ! all models available will be included in ensemble.modeling")
    chosen.models <- modeling.output@models.computed
  } else{
    chosen.models.check <- chosen.models %in% modeling.output@models.computed
    if(sum(!chosen.models.check) > 0){
      stop(paste("Some selected models do not exist: ", toString(chosen.models[!chosen.models.check]),
                 "\nPlease choose models in computed models ( ",
                 toString(modeling.output@models.computed), " )",sep=""))
    }
  }

  # 3. eval.metric checking
  if(!is.null(eval.metric)){
    if(!is.character(eval.metric)){
      stop("eval.metric must be a character vector or NULL")
    }
    if('all' %in% eval.metric){
      eval.metric <- dimnames(get_evaluations(modeling.output))[[1]]
    }
    eval.metric.check <- eval.metric %in% dimnames(get_evaluations(modeling.output))[[1]]
    if(sum(!eval.metric.check) > 0){
      stop(paste("Some selected evaluation metrics are not available: ", toString(eval.metric[!eval.metric.check]),
                 "\nPlease choose some in those computed yet ( ",
                 toString(dimnames(get_evaluations(modeling.output))[[1]]), " )",sep=""))
    }
  }

  # 4. eval.metric.quality.threshold
  if(!is.null(eval.metric)){
    if(!is.null(eval.metric.quality.threshold)){
      if(!is.numeric(eval.metric.quality.threshold)){
        stop("eval.metric.quality.threshold must be NULL or a numeric vector")
      }
      if(length(eval.metric) != length(eval.metric.quality.threshold)){
        stop("you must specify as many eval.metric.quality.threshold as eval.metric (if you specify some)")
      }
      cat("\n   > Evaluation & Weighting methods summary :\n")
      cat(paste(eval.metric, eval.metric.quality.threshold,  sep = " over ", collapse = "\n      "), fill=TRUE, labels = "     ")
    } else{
      cat("\n   ! No eval.metric.quality.threshold -> All models will be kept for Ensemble Modeling")
      eval.metric.quality.threshold <- rep(0, length(eval.metric))
    }
  }

  # 4b. model.eval.meth checking
  models.eval.meth <- unique(models.eval.meth)

  if(sum(models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS','KAPPA','ACCURACY','BIAS',
                                 'POD','PODFD','CSI','ETS','HK','ROC')) != length(models.eval.meth)){
    stop(paste(models.eval.meth[which( (models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS',
                                                                'KAPPA','ACCURACY','BIAS', 'POD',
                                                                'PODFD','CSI', 'ETS','HK','ROC'))
                                       == FALSE) ]," is not an availabe model evaluation metric !",sep=""))
  }

  # 5. check selected EM algo
  if( !is.logical(prob.mean) | !is.logical(prob.cv) | !is.logical(prob.ci) | !is.logical(prob.median) |
      !is.logical(committee.averaging) | !is.logical(prob.mean.weight) ){
    stop("prob.mean, prob.cv, prob.ci, prob.median, committee.averaging and prob.mean.weight arguments must be logical")
  }
  if(is.null(eval.metric)){
    if(committee.averaging | prob.mean.weight){
      stop("You must choose eval.metric if you want to compute Committee Averaging or Probability weighting mean algorithms")
    }
  }

  # 6. alpha for Confident interval
  if(prob.ci){
    if(!is.numeric(prob.ci.alpha)){
      stop("prob.ci.alpha must be numeric")
    }
    if(prob.ci.alpha <= 0 | prob.ci.alpha>= 0.5){
      stop("prob.ci.alpha must be a numeric between 0 and 0.5")
    }
  }

  # 7. decay checking
  if(prob.mean.weight){
    test.prob.mean.weight.decay <- TRUE
    ## check compatibility of prob.mean.weight.decay class
    if(!is.numeric(prob.mean.weight.decay) & !is.character(prob.mean.weight.decay) & !is.function(prob.mean.weight.decay)){
      test.prob.mean.weight.decay <- FALSE
    } else if(is.numeric(prob.mean.weight.decay)){ ## check numeric prob.mean.weight.decay
      if(prob.mean.weight.decay < 0){
        test.prob.mean.weight.decay <- FALSE
      }
    } else if(is.character(prob.mean.weight.decay)){ ## check character prob.mean.weight.decay
      if(prob.mean.weight.decay != 'proportional'){
        test.prob.mean.weight.decay <- FALSE
      }
    }

    if(!test.prob.mean.weight.decay){
      stop("'prob.mean.weight.decay' should be either 'proportional', a numeric value > 0 or a function")
    }
  }

  if(is.null(eval.metric)){
    eval.metric <- 'none'
  }

  # 8. by arg checking
  available.em.by <- c('PA_dataset', 'algo', 'all', 'PA_dataset+repet', 'PA_dataset+algo')
  if(!(em.by %in% available.em.by) ){
    stop("Invalid 'em.by' argument given. It must be one of : 'PA_dataset', 'algo', 'all', 'PA_dataset+repet' or 'PA_dataset+algo'")
  }


  return( list( modeling.output = modeling.output,
                chosen.models = chosen.models,
                eval.metric = eval.metric,
                eval.metric.quality.threshold = eval.metric.quality.threshold,
                models.eval.meth = models.eval.meth,
                prob.mean = prob.mean,
                prob.cv = prob.cv,
                prob.ci = prob.ci,
                prob.ci.alpha = prob.ci.alpha,
                prob.median = prob.median,
                committee.averaging = committee.averaging,
                prob.mean.weight = prob.mean.weight,
                prob.mean.weight.decay  = prob.mean.weight.decay,
                em.by = em.by))

}


# =-=-=-=-=-=-=-=- em.models.assembling function -=-=-=-=-=-=-=- #
.em.models.assembling <- function(chosen.models, em.by){
  assembl.list = list()

  if(em.by == 'PA_dataset'){
    for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
#       assembl.list[[paste(dat,"_AllRun", sep="")]] <- chosen.models[grep(paste("_",dat,"_",sep=""), chosen.models)]
      assembl.list[[paste("mergedAlgo_mergedRun_", dat, sep="")]] <- chosen.models[grep(paste("_",dat,"_",sep=""), chosen.models)]
    }
    return(assembl.list)
  }

  if(em.by == 'algo'){
    for(algo in .extractModelNamesInfo(chosen.models, info='models')){
#       assembl.list[[paste(algo,"_AllRun", sep="")]] <- chosen.models[grep(paste("_",algo,sep=""), chosen.models)]
      assembl.list[[paste(algo,"_mergedRun_mergedData", sep="")]] <- chosen.models[grep(paste("_",algo,sep=""), chosen.models)]
    }
    return(assembl.list)
  }

  if(em.by == 'all'){
#     assembl.list[["TotalConsensus"]] <- chosen.models
    assembl.list[[paste("mergedAlgo_mergedRun_mergedData", sep="")]] <- chosen.models
    return(assembl.list)
  }

  if(em.by == 'PA_dataset+repet'){
    for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
      for(repet in .extractModelNamesInfo(chosen.models, info='run.eval')){
        mod.tmp <- intersect(x=grep(paste("_",dat,"_",sep=""), chosen.models),
                             y=grep(paste("_",repet,"_",sep=""), chosen.models))
        if(length(mod.tmp)){
#           assembl.list[[paste(dat,"_",repet,'_AllAlgos', sep="")]] <- chosen.models[mod.tmp]
          assembl.list[[paste("mergedAlgo_",repet,"_",dat, sep="")]] <- chosen.models[mod.tmp]
        }
      }
    }
    return(assembl.list)
  }

  if(em.by == 'PA_dataset+algo'){
    for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
      for(algo in .extractModelNamesInfo(chosen.models, info='models')){
        mod.tmp <- intersect(x=grep(paste("_",dat,"_",sep=""), chosen.models),
                             y=grep(paste("_",algo,sep=""), chosen.models))
        if(length(mod.tmp)){
#           assembl.list[[paste(dat,"_AllRepet_",algo, sep="")]] <- chosen.models[mod.tmp]
          assembl.list[[paste(algo,"_mergedRun_", dat, sep="")]] <- chosen.models[mod.tmp]
        }
      }
    }
    return(assembl.list)
  }

}
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



.get_needed_predictions <- function(modeling.output, em.by, models.kept, eval.metric, eval.metric.quality.threshold){
  out <- list(predictions = NULL,
              models.kept = NULL,
              models.kept.scores = NULL)
  for(eval.m in eval.metric){
    if( eval.m != 'none'){
      models.kept.scores <- unlist(lapply(models.kept, function(x){
        mod <- tail(unlist(strsplit(x,"_")), 3)[3]
        run <- tail(unlist(strsplit(x,"_")), 3)[2]
        dat <- tail(unlist(strsplit(x,"_")), 3)[1]
        # select evaluations scores obtained for Evaluation Data if exists or CV if not
        if(modeling.output@has.evaluation.data){
          return(get_evaluations(modeling.output)[eval.m, "Evaluating.data", mod, run, dat])
        } else{
          return(get_evaluations(modeling.output)[eval.m, "Testing.data", mod, run, dat])
        }

      }))
      ## set NA to -1
      if(!is.null(models.kept.scores)){
        models.kept.scores[is.na(models.kept.scores)] <- -1
      }
      out$models.kept[[eval.m]] <- models.kept[models.kept.scores > eval.metric.quality.threshold[which(eval.metric == eval.m)]]
      out$models.kept.scores[[eval.m]] <- models.kept.scores[models.kept.scores > eval.metric.quality.threshold[which(eval.metric == eval.m)]]
    } else{
      out$models.kept[[eval.m]] <- models.kept
    }
  }

  models.kept.union <- unique(unlist(out$models.kept))

  if(length(models.kept.union) ){
#     if(modeling.output@has.evaluation.data){
#       out$predictions <- as.data.frame(get_predictionsEval(modeling.output, as.data.frame = TRUE)[,models.kept.union, drop=F])
#     } else{
      ## load prediction on each PA dataset
      if(em.by %in% c("PA_dataset",'PA_dataset+algo','PA_dataset+repet')){
        out$predictions <- as.data.frame(get_predictions(modeling.output, as.data.frame = TRUE)[,models.kept.union, drop=F])
      } else{ ## redo prediction on full data.set
        cat("\n   ! Models projections for whole zonation required...")
        temp_name <- paste('tmp_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep="")
        out$predictions <- BIOMOD_Projection(modeling.output = modeling.output,
                                             new.env = get_formal_data(modeling.output)@data.env.var,
                                             proj.name = temp_name,
                                             xy.new.env = get_formal_data(modeling.output)@coord,
                                             selected.models = models.kept.union,
                                             compress = TRUE,
                                             build.clamping.mask = F,
                                             do.stack=T, silent = T)@proj@val

        # transform array into data.frame
        out$predictions <- as.data.frame(out$predictions)
        names(out$predictions) <- unlist(lapply(strsplit(names(out$predictions),".", fixed=TRUE),
                                         function(x){
                                           x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                           data.set.id <- x.rev[1]
                                           cross.valid.id <- x.rev[2]
                                           algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                           model.id <- paste(modeling.output@sp.name,
                                                             data.set.id,
                                                             cross.valid.id,
                                                             algo.id, sep="_")
                                           return(model.id)
                                         }))
        # keep only wanted columns
        out$predictions <- out$predictions[,models.kept.union, drop=F]
        unlink(file.path(modeling.output@sp.name,paste("proj_", temp_name, sep="") ),recursive = TRUE, force = TRUE)
        cat("\n")
      }
#     }
    return(out)
  } else {
    cat("\n   ! No models kept due to threshold filtering... Ensemble Modeling was skipped!")
    return(NULL)
  }
}




####################################################################################################
# BIOMOD_FormatingData
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   puting input data in the right format and doing Pseudo Absences selection if desire

# INPUT :
#   resp.var <- Response Variable (monospecific) as vector, sp.point.data.frame or rasterLayer
#               code for vector and sp.objects : 1=pres, 0=true_abs, NA=no_info
#               code for vector and sp.objects : 1=pres, 0=true_abs, -1=no_info
#   expl.var <- Explanatory Variable as matrix, data.frame, sp.point.data.frame or rasterStack
#   resp.xy <- coordiantes of reponse points (2 column matrix)
#   resp.name <- name of considered specie
#   eval.resp.var <- independent response variable for models evaluations
#   eval.expl.var <- independent explanatory variable for models evaluations
#   eval.resp.xy <- independent response variable coordinates variable for models evaluations
#   PA.nb.rep <- Nb of Pseudo Absences Run to compute
#   PA.nb.absences <- Nb of Absences selected (true absences are counted in)
#   PA.strategy <- Pseudo Absences strategy
#   PA.dist.min <- If strategy is 'disk' : Pseudo Absences minimum distance between pres and selected absences (in metters if explanatory is georeferenced or in resp.xy units in all other cases)
#   PA.dist.man <- If strategy is 'disk' : Pseudo Absences maximum distance between pres and selected absences (in metters if explanatory is georeferenced or in resp.xy units in all other cases)
#
#   PA.sre.quant <- If strategy is 'sre' : the quantile use for sre calculation
#   PA.table <- If strategy is 'user.defined' : a boolean data.frame indiacating which points of resp.var should be sonsidered in each PA run.
#   na.rm <- if True na are automatically removed

# OUTPUT :
#   a BIOMOD.formated.data object that will be given to BIOMOD_Modeling function

####################################################################################################

'BIOMOD_FormatingData' <- function(resp.var,
                                   expl.var,
                                   resp.xy = NULL,
                                   resp.name = NULL,
                                   eval.resp.var = NULL,
                                   eval.expl.var = NULL,
                                   eval.resp.xy = NULL,
                                   PA.nb.rep = 0,
                                   PA.nb.absences = 1000,
                                   PA.strategy = 'random',
                                   PA.dist.min = 0,
                                   PA.dist.max = NULL,
                                   PA.sre.quant = 0.025,
                                   PA.table = NULL,
                                   na.rm = TRUE){
  .bmCat(paste(resp.name, " Data Formating", sep=""))

  # 1 check args
  args <- .BIOMOD_FormatingData.check.args(resp.var,
                                           expl.var,
                                           resp.xy,
                                           resp.name,
                                           eval.resp.var,
                                           eval.expl.var,
                                           eval.resp.xy,
                                           PA.nb.rep,
                                           PA.nb.absences,
                                           PA.strategy,
                                           PA.dist.min,
                                           PA.dist.max,
                                           PA.sre.quant,
                                           PA.table)

  resp.var <- args$resp.var
  expl.var <- args$expl.var
  resp.xy <- args$resp.xy
  resp.name <- args$resp.name
  eval.resp.var <- args$eval.resp.var
  eval.expl.var <- args$eval.expl.var
  eval.resp.xy <- args$eval.resp.xy
  PA.nb.rep <- args$PA.nb.rep
  PA.nb.absences <- args$PA.nb.absences
  PA.strategy <- args$PA.strategy
  PA.dist.min <- args$PA.dist.min
  PA.dist.max <- args$PA.dist.max
  PA.sre.quant <- args$PA.sre.quant
  PA.table <- args$PA.table

  rm('args')
  gc()

  out <- NULL

  if(PA.strategy == 'none'){ # no Pseudo Absences
    out <- BIOMOD.formated.data(sp=resp.var,
                                xy=resp.xy,
                                env=expl.var,
                                sp.name=resp.name,
                                eval.sp=eval.resp.var,
                                eval.env=eval.expl.var,
                                eval.xy=eval.resp.xy,
                                na.rm=na.rm)
  } else{ # Automatic Pseudo Absences Selection
    out <- BIOMOD.formated.data.PA(sp=resp.var, xy=resp.xy, env=expl.var, sp.name=resp.name,
                                   eval.sp=eval.resp.var, eval.env=eval.expl.var, eval.xy=eval.resp.xy,
                                   PA.NbRep=PA.nb.rep, PA.strategy=PA.strategy,
                                   PA.nb.absences = PA.nb.absences, PA.dist.min = PA.dist.min,
                                   PA.dist.max = PA.dist.max, PA.sre.quant = PA.sre.quant, PA.table=PA.table,
                                   na.rm=na.rm)
  }


  .bmCat("Done")
  return(out)
}

.BIOMOD_FormatingData.check.args <- function(resp.var,
                                             expl.var,
                                             resp.xy,
                                             resp.name,
                                             eval.resp.var,
                                             eval.expl.var,
                                             eval.resp.xy,
                                             PA.nb.rep,
                                             PA.nb.absences,
                                             PA.strategy,
                                             PA.dist.min,
                                             PA.dist.max,
                                             PA.sre.quant,
                                             PA.table){

  # 0. names checking



  ### check resp.name is available
  if(grepl('_',resp.name) | grepl(' ',resp.name)){
    resp.name <- paste(unlist(strsplit(resp.name,'_')),collapse='.')
    resp.name <- paste(unlist(strsplit(resp.name,' ')),collapse='.')

    cat('\n Response variable name was converted into', resp.name)
  }

  ### check resp.name is available
  ### Not done because no imporance

  # 1. Checking input params class
  available.types <- c( 'numeric', 'data.frame', 'matrix',
                        'RasterLayer', 'RasterStack',
                        'SpatialPointsDataFrame', 'SpatialPoints')
  ###### resp.var
  if(!inherits(resp.var, available.types)){
    stop( paste("Response variable must be one of ", toString(available.types), sep="") )
  }

  ### response var raster object not supported yet
  if(inherits(resp.var, 'Raster')){
    stop("Raster response variable not supported yet ! \nPlease extract your Presences and your absences by yourself")
    #### TO DO ####
    ## extract the 0 and 1 in sp format
  }

  ###### expl.var
  if(!inherits(expl.var, setdiff(available.types, 'SpatialPoints'))){
    stop( paste("Explanatory variable must be one of ", toString(available.types), sep="") )
  }


  ###### resp.xy
  if(inherits(resp.var,'SpatialPoints') ){
    if(!is.null(resp.xy)){
      cat("\n      ! XY coordinates of response variable will be ignored because spatial response object is given.")
    }
    resp.xy <- data.matrix(sp::coordinates(resp.var))
    if(inherits(resp.var, 'SpatialPointsDataFrame')){
      resp.var <- resp.var@data
    } else{
      cat("\n      ! Response variable is considered as only presences... Is it really what you want?")
      resp.var <- rep(1,nrow(resp.xy))
    }

  }


  ### transforming into numeric if data.frame or matrix
  if(is.matrix(resp.var) | is.data.frame(resp.var)){
    if(ncol(resp.var) > 1){
      stop("You must give a monospecific response variable (1D object)")
    } else{
      resp.var <- as.numeric(resp.var[,1])
    }
  }

  if(is.matrix(expl.var) | is.numeric(expl.var) ){
    expl.var <- as.data.frame(expl.var)
  }

  if(inherits(expl.var, 'Raster')){
    expl.var <- raster::stack(expl.var, RAT=FALSE)
  }

  if(inherits(expl.var, 'SpatialPoints')){
    expl.var <- as.data.frame(expl.var@data)
  }

  ### check of xy coordinates validity
  if(!is.null(resp.xy)){
    if(ncol(resp.xy)!=2){
      stop("if given, resp.xy must be a 2 column matrix or data.frame")
    }
    if(nrow(resp.xy) != length(resp.var)){
      stop("Response variable and its coordinates don't match")
    }
    resp.xy <- as.data.frame(resp.xy)
  }

  ### convert response var into binary
  resp.var[which(resp.var>0)] <- 1
  resp.var[which(resp.var<=0)] <- 0

  #### At this point :
  ####  - resp.var is a numeric
  ####  - resp.xy is NULL or a data.frame
  ####  - expl.var is a data.frame or a RasterStack
  ####  - sp.name is a character

  ### check resp and expl var compatibility
  if(is.data.frame(expl.var)){
    if(nrow(expl.var) != length(resp.var)){
      stop("If explanatory variable is not a raster then dimensions of response variable and explanatory variable must match!")
    }
  }

  ### PA strategy
#   if(!is.null(PA.strategy)){ # force PA.nb.rep to be positive if PA.strategy is defined
#     PA.nb.rep = max(c(PA.nb.rep,1))
#   }

  if(is.null(PA.table) & PA.nb.rep < 1){
    cat("\n> No pseudo absences selection !")
    PA.strategy <- "none"
    PA.nb.rep <- 0
  }

  if(is.null(PA.strategy) &  PA.nb.rep > 0){
    cat("\n> Pseudo absences will be selected randomly !")
    PA.strategy <- "random"
  }


  if( !is.null(PA.table)){
    cat("\n> Pseudo absences used will be user defined ones !")
    PA.strategy <- "user.defined"
    PA.nb.rep <- 0
  }

  if(PA.strategy == "user.defined"){
    if(! (is.matrix(PA.table) | is.data.frame(PA.table)))
      stop("\n PA.table must be a matrix or a data.frame")

    if(nrow(PA.table) != length(resp.var))
      stop("\n PA.table must have as many row than the number
           of observation of your response variable")

    #PA.table <- as.data.frame(sapply(PA.table,simplify=FALSE,as.logical))
    colnames(PA.table) <- paste("PA",1:ncol(PA.table),sep="")

  }

  # 2. eval.resp.var.checking

  if(!is.null(eval.resp.var)){
    # do the same test than previous one
    ###### eval.resp.var
    if(!(class(eval.resp.var) %in% available.types)){
      stop( paste("Response variable must be one of ", toString(available.types), sep="") )
    }

    ### response var raster object not supported yet
    if(inherits(eval.resp.var, 'Raster')){
      stop("Raster response variable not supported yet ! \nPlease extract your Presences and your absences by yourself")
      #### TO DO ####
      ## extract the 0 and 1 in sp format
    }

    ###### expl.var
    if(!is.null(eval.expl.var)){
      if(!(class(eval.expl.var) %in% available.types[-which(available.types == 'SpatialPoints')])){
        stop( paste("Explanatory variable must be one of ", toString(available.types), sep="") )
      }
    } else{
      if(!(inherits(expl.var, 'Raster'))){
        stop("If explanatory variable is not a raster and you want to consider evaluation response variable, you have to give evaluation explanatory variables")
      }
    }

    ###### resp.xy
    if(inherits(eval.resp.var,'SpatialPoints') ){
      if(!is.null(eval.resp.xy)){
        cat("\n      ! XY coordinates of response variable will be ignored because spatial response object is given.")
      }
      eval.resp.xy <- data.matrix(sp::coordinates(eval.resp.var))
      if(class(eval.resp.var) == 'SpatialPointsDataFrame'){
        eval.resp.var <- eval.resp.var@data
      } else{
        cat("\n      ! Response variable is considered as only presences... Is it really what you want?")
        eval.resp.var <- rep(1,nrow(eval.resp.xy))
      }

    }


    ### transforming into numeric if data.frame or matrix
    if(is.matrix(eval.resp.var) | is.data.frame(eval.resp.var)){
      if(ncol(eval.resp.var) > 1){
        stop("You must give a monospecific response variable (1D object)")
      } else{
        eval.resp.var <- as.numeric(eval.resp.var[,1])
      }
    }

    if(is.matrix(eval.expl.var) | is.numeric(eval.expl.var) ){
      eval.expl.var <- as.data.frame(eval.expl.var)
    }

    if(inherits(eval.expl.var, 'Raster')){
      eval.expl.var <- raster::stack(eval.expl.var)
    }

    if(inherits(eval.expl.var, 'SpatialPoints')){
      eval.expl.var <- as.data.frame(eval.expl.var@data)
    }

    ### check of xy coordinates validity
    if(!is.null(eval.resp.xy)){
      if(ncol(eval.resp.xy)!=2){
        stop("if given, resp.xy must be a 2 column matrix or data.frame")
      }
      if(nrow(eval.resp.xy) != length(eval.resp.var)){
        stop("Response variable and its coordinates don't match")
      }
      eval.resp.xy <- as.data.frame(eval.resp.xy)
    }

    if(is.data.frame(eval.expl.var)){
      if(nrow(eval.expl.var) != length(eval.resp.var)){
        stop("If explanatory variable is not a raster then dimensions of response variable and explanatory variable must match!")
      }
    }

    ### remove NAs from evaluation data
    if( sum(is.na(eval.resp.var)) > 0 ){
      cat("\n      ! NAs have been automatically removed from Evaluation data")
      if(!is.null(eval.resp.xy)){
        eval.resp.xy <- eval.resp.xy[-which(is.na(eval.resp.var)),]
      }
      eval.resp.var <- na.omit(eval.resp.var)
    }

    ### convert response var into binary
    eval.resp.var[which(eval.resp.var>0)] <- 1
    eval.resp.var[which(eval.resp.var<=0)] <- 0

    ### check there are both presences and absences in evaluation dataset
    if( sum(eval.resp.var == 1) < 1 | sum(eval.resp.var == 0) < 1){
      stop("Evaluation response data must have both presences and absences")
    }

  } else {
    cat("\n      ! No data has been set aside for modeling evaluation")
    eval.expl.var <- eval.resp.xy <- NULL
  }

  ### PA arguments are not checked here because it will be done later... (may be will do it here later)

  return(list( resp.var = resp.var,
               expl.var = expl.var,
               resp.xy = resp.xy,
               resp.name = resp.name,
               eval.resp.var = eval.resp.var,
               eval.expl.var = eval.expl.var,
               eval.resp.xy = eval.resp.xy,
               PA.nb.rep = PA.nb.rep,
               PA.nb.absences = PA.nb.absences,
               PA.strategy = PA.strategy,
               PA.dist.min = PA.dist.min,
               PA.dist.max = PA.dist.max,
               PA.sre.quant = PA.sre.quant,
               PA.table = PA.table))

}
BIOMOD_LoadModels <- function(bm.out, ... ){

  ####
  # ... can be models, run.eval, data.set to make a models subselection
  add.args <- list(...)
  args <- .BIOMOD_LoadModels.check.args(bm.out, add.args)
  add.args <- args$add.args
  rm(args)

  models.to.load <- get_built_models(bm.out)

  envir <- parent.frame()


  if(!is.null(add.args$full.name)){
    models.to.load <- add.args$full.name
  } else{ ## make a subselection
    if(!is.null(add.args$models)){
      model.to.load.tmp <- c()
      for(mod in add.args$models){
        if(sum(grepl(mod,  models.to.load)) > 0){
          model.to.load.tmp <- c(model.to.load.tmp, grep(mod, models.to.load, value=TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }

    if(!is.null(add.args$run.eval)){
      model.to.load.tmp <- c()
      for(re in add.args$run.eval){
        if(sum(grepl(re,  models.to.load)) > 0){
          model.to.load.tmp <- c(model.to.load.tmp, grep(re, models.to.load, value=TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }

    if(!is.null(add.args$data.set)){
      model.to.load.tmp <- c()
      for(ds in add.args$data.set){
        if(sum(grepl(ds,  models.to.load)) > 0){
          model.to.load.tmp <- c(model.to.load.tmp, grep(ds, models.to.load, value=TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }

  }

  if(length(models.to.load) == 0){
    cat("\n   ! No models computed matched, No models loaded !")
    return(NULL)
  }

  if(!is.null(add.args$as) & length(models.to.load)==1){
    assign(x=add.args$as,
           value=get(load(file=file.path(bm.out@sp.name, "models", bm.out@modeling.id, models.to.load))),
           envir=envir)
    invisible(TRUE)
  } else {
    for(mtl in models.to.load){
      load(file=file.path(bm.out@sp.name, "models", bm.out@modeling.id, mtl),envir=envir)
    }

    return(models.to.load)
  }


}

.BIOMOD_LoadModels.check.args <- function(bm.out, add.args){
  if(! (inherits(bm.out, 'BIOMOD.models.out') | inherits(bm.out, 'BIOMOD.EnsembleModeling.out'))){
    stop("bm.out arg must be a BIOMOD.models.out or a BIOMOD.EnsembleModeling.out object")
  }

  available.args <- c("models", "run.eval", "data.set", "path", "as", "full.name")
  given.args <- names(add.args)

  ## is all additional args are known ?
  if(sum(given.args %in% available.args) != length(given.args)){
    cat("\n   !", toString( given.args[which(!(given.args %in% available.args))] ), "arguments are unknown. Please refer to function help file to see available ones.", fill=.Options$width)
    ## remove unknown args
    for(ga in given.args[which(!(given.args %in% available.args))]){
      add.args[[ga]] <- NULL
    }
  }

  ## get all available model names
  avail_models <- get_built_models(bm.out)

  ## check additional args values
  ### models names
  if(!is.null(add.args$models)){
    if(sum(add.args$models %in% .extractModelNamesInfo(model.names=avail_models, info='models') ) != length(add.args$models) ){
      stop(paste("models argument must be one of ", toString(.extractModelNamesInfo(model.names=avail_models, info='models')), sep="") )
    } else{
      add.args$models = paste("_", add.args$models, sep="")
    }
  }

  ### run.eval names
  if(!is.null(add.args$run.eval)){
    if(sum(add.args$run.eval %in% .extractModelNamesInfo(model.names=avail_models, info='run.eval') != length(add.args$run.eval)) ){
      stop(paste("run.eval argument must be one of ", toString(.extractModelNamesInfo(model.names=avail_models), info='run.eval'), sep="") )
    } else{
      add.args$run.eval = paste("_", add.args$run.eval, "_", sep="")
    }
  }

  ### data.set names
  if(!is.null(add.args$data.set)){
    if(sum(add.args$data.set %in% .extractModelNamesInfo(model.names=avail_models, info='data.set') != length(add.args$data.set)) ){
      stop(paste("data.set argument must be one of ", toString(.extractModelNamesInfo(model.names=avail_models), info='data.set'), sep="") )
    } else{
      add.args$data.set = paste("_", add.args$data.set, "_", sep="")
    }
  }

  ### path to sim
  if(!is.null(add.args$path)){
    if(!(bm.out@sp.name %in% list.dirs(path = add.args$path))){
      stop("invalid path given")
    }
  } else{
    add.args$path = "."
  }

  ### full.name
  if(!is.null(add.args$full.name)){
    if(sum(!(add.args$full.name %in% avail_models)) > 0){
      stop("full.name arg must be one of : ", toString(avail_models))
    }
  }


  return(list(add.args = add.args))

}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

'.extractModelNamesInfo' <- function(model.names, info = 'species'){
  if(!is.character(model.names)){
    stop("model.names must be a character vector")
  }
  if(!is.character(info) | length(info) != 1 | !(info %in% c('species', 'data.set', 'models', 'run.eval')) ){
    stop("info must be 'species', 'data.set', 'models' or 'run.eval'")
  }

  info.tmp <- as.data.frame(strsplit(model.names, "_"))

  return( switch(info,
                 species = paste(unique(unlist(info.tmp[-c(nrow(info.tmp), nrow(info.tmp)-1, nrow(info.tmp)-2),])), collapse="_"),
                 data.set = paste(unique(unlist(info.tmp[(nrow(info.tmp)-2),]))),
                 run.eval = paste(unique(unlist(info.tmp[(nrow(info.tmp)-1),]))),
                 models = paste(unique(unlist(info.tmp[(nrow(info.tmp)),]))) ) )

}
##' @name BIOMOD_Modeling
##' @aliases BIOMOD_Modeling
##' @title Run a range of species distribution models
##' @description
##' This function allows to calibrate and evaluate a range of
##' species distribution models techniques run over a given
##' species. Calibrations are made on the whole sample or a
##' random subpart. The predictive power of the different models
##' is estimated using a range of evaluation metrics.
##' 
##' @param data \code{BIOMOD.formated.data} object returned by
##'   \code{\link[biomod2]{BIOMOD_FormatingData}}
##' @param models character, models to be computed names. To be 
##'   chosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE',
##'   'FDA', 'MARS', 'RF', 'MAXENT.Phillips', 'MAXENT.Phillips.2'
##' @param models.options \code{BIOMOD.models.options} object
##'   returned by \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' @param NbRunEval integer, number of Evaluation run.
##' @param DataSplit numeric, \% of data used to calibrate the
##'   models, the remaining part will be used for testing
##' @param Yweights numeric, vector of weights (one per 
##'   observation)
##' @param Prevalence either \code{NULL} (default) or a 0-1
##'   numeric used to build 'weighted response weights'
##' @param VarImport Number of permutation to estimate variable
##'   importance
##' @param models.eval.meth vector of names of evaluation metric
##'   among 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY',
##'   'BIAS', 'POD', 'CSI' and 'ETS'
##' @param SaveObj keep all results and outputs on hard drive or
##'   not (NOTE: strongly recommended)
##' @param rescal.all.models if true, all model prediction will
##'   be scaled with a binomial GLM
##' @param do.full.models if true, models calibrated and
##'   evaluated with the whole dataset are done
##' @param modeling.id character, the ID (=name) of modeling
##'   procedure. A random number by default.
##' @param \ldots further arguments :
##' 
##'  - \code{DataSplitTable} : a \code{matrix}, \code{data.frame}
##'    or a 3D \code{array} filled with \code{TRUE/FALSE} to
##'    specify which part of data must be used for models
##'    calibration (\code{TRUE}) and for models validation
##'    (\code{FALSE}). Each column corresponds to a 'RUN'. If
##'    filled, args \code{NbRunEval}, \code{DataSplit} and
##'    \code{do.full.models} will be ignored.
##'    
##' @details 
##' 
##' 1. \bold{data}
##' .. If you have decide to add pseudo absences to your
##' original dataset (see 
##' \code{\link[biomod2]{BIOMOD_FormatingData}}), 
##' NbPseudoAbsences * \code{NbRunEval + 1} models will be
##' created.
##' 
##' 2. \bold{models}
##' .. The set of models to be calibrated on the data. 10
##' modeling techniques are currently available:
##' 
##' .. - GLM : Generalized Linear Model 
##' (\code{\link[stats]{glm}})
##' 
##' .. - GAM : Generalized Additive Model (\code{\link[gam]{gam}},
##' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}, see 
##' \code{\link[biomod2]{BIOMOD_ModelingOptions} for details on 
##' algorithm selection})
##' 
##' .. - GBM : Generalized Boosting Model or usually called Boosted
##' Regression Trees (\code{\link[gbm]{gbm}})
##' 
##' .. - CTA: Classification Tree Analysis (\code{\link[rpart]{rpart}})
##' 
##' .. - ANN: Artificial Neural Network (\code{\link[nnet]{nnet}})
##' 
##' .. - SRE: Surface Range Envelop or usually called BIOCLIM
##' 
##' .. - FDA: Flexible Discriminant Analysis (\code{\link[mda]{fda}})
##' 
##' .. - MARS: Multiple Adaptive Regression Splines 
##' (\code{\link[earth]{earth}})
##' 
##' .. - RF: Random Forest (\code{\link[randomForest]{randomForest}})
##' 
##' .. - MAXENT.Phillips: Maximum Entropy (
##' \url{https://biodiversityinformatics.amnh.org/open_source/maxent})
##' 
##' .. - MAXENT.Phillips.2: Maximum Entropy 
##' (\code{\link[maxnet]{maxnet}})
##' 
##' 3. \bold{NbRunEval & DataSplit}
##' .. As already explained in the \code{\link{BIOMOD_FormatingData}}
##' help file, the common trend is to split the original dataset into 
##' two subsets, one to calibrate the models, and another one to evaluate
##' them. Here we provide the possibility to repeat this process
##' (calibration and evaluation) N times (\code{NbRunEval} times). 
##' The proportion of data kept for calibration is determined by the
##' \code{DataSplit} argument (100\% - \code{DataSplit} will be used to
##' evaluate the model). This sort of cross-validation allows to have a
##' quite robust test of the models when independent data are not
##' available. Each technique will also be calibrated on the complete
##' original data. All the models produced by BIOMOD and their related
##' informations are saved on the hard drive.
##' 
##' 4. \bold{Yweights & Prevalence}
##' .. Allows to give more or less weight to some particular 
##' observations. If these arguments is kept to NULL 
##' (\code{Yweights = NULL}, \code{Prevalence = NULL}), each 
##' observation (presence or absence) has the same weight (independent 
##' of the number of presences and absences). If \code{Prevalence = 0.5} 
##' absences will be weighted equally to the presences (i.e. the 
##' weighted sum of presence equals the weighted sum of absences). If
##' prevalence is set below or above 0.5 absences or presences are given
##' more weight, respectively.
##' .. In the particular case that pseudo-absence data have been
##' generated \code{BIOMOD_FormatingData} (\code{PA.nb.rep > 0}), weights
##' are by default (\code{Prevalence = NULL}) calculated such that
##' prevalence is 0.5, meaning that the presences will have the same
##' importance as the absences in the calibration process of the models.
##' Automatically created \code{Yweights} will be composed of integers to
##' prevent different modeling issues.
##' .. Note that the \code{Prevalence} argument will always be ignored if
##' \code{Yweights} are defined.
##' 
##' 5. \bold{models.eval.meth}
##' .. The available evaluations methods are :
##' 
##' .. - \code{ROC} : Relative Operating Characteristic
##' .. - \code{KAPPA} : Cohen's Kappa (Heidke skill score)
##' .. - \code{TSS} : True kill statistic (Hanssen and Kuipers 
##' discriminant, Peirce's skill score)
##' .. - \code{FAR} : False alarm ratio
##' .. - \code{SR} : Success ratio
##' .. - \code{ACCURANCY} : Accuracy (fraction correct)
##' .. - \code{BIAS} : Bias score (frequency bias)
##' .. - \code{POD} : Probability of detection (hit rate)
##' .. - \code{CSI} : Critical success index (threat score)
##' .. - \code{ETS} : Equitable threat score (Gilbert skill score)
##' 
##' Some of them are scaled to have all an optimum at 1. You can choose
##' one of more (vector) evaluation metric. By Default, only 'KAPPA',
##' 'TSS' and 'ROC' evaluation are done. Please refer to the CAWRC
##' website (\url{http://www.cawcr.gov.au/projects/verification/##'Methods_for_dichotomous_forecasts}) 
##' to get detailed description of each metric.
##' 
##' 6. \bold{SaveObj}
##' If this argument is set to False, it may prevent the evaluation of
##' the \sQuote{ensemble modeled} models in further steps. We strongly
##' recommend to always keep this argument \code{TRUE} even it asks for
##' free space onto the hard drive.
##'
##' 7. \bold{rescal.all.models}
##' \bold{This parameter is quite experimental and we advise not to use
##' it. It should lead to reduction in projection scale amplitude}
##' Some categorical models have to be scaled in every case (
##' \sQuote{FDA}, \sQuote{ANN}). But It may be interesting to scale all
##' model computed to ensure that they will produced comparable 
##' predictions (0-1000 ladder). That's particularly useful when you 
##' do some ensemble forecasting to remove the scale prediction effect
##' (the more extended projections are, the more they influence ensemble
##' forecasting results).
##' 
##' 8. \bold{do.full.models}
##' Building models with all information available may be useful in some
##' particular cases (i.e. rare species with few presences points). The
##' main drawback of this method is that, if you don't give separated
##' data for models evaluation, your models will be evaluated with the
##' same data that the ones used for calibration. That will lead to 
##' over-optimistic evaluation scores. Be careful with this '_Full'
##' models interpretation.
##' 
##' @return
##' A BIOMOD.models.out object
##' See \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} 
##' for details.
##' Additional objects are stored out of R in two different directories
##' for memory storage purposes. They are created by the function
##' directly on the root of your working directory set in R ("models"
##' directory). This one contains each calibrated model for each
##' repetition and pseudo-absence run. A hidden folder 
##' \code{.DATA_BIOMOD} contains some files (predictions, original
##' dataset copy, pseudo absences chosen...) used by other functions like
##' \code{\link[biomod2]{BIOMOD_Projection}} or 
##' \code{\link[biomod2]{BIOMOD_EnsembleModeling}}.
##' 
##' The models are currently stored as objects to be read exclusively in
##' R. To load them back (the same stands for all objects stored on the
##' hard disk) use the \code{\link{load}} function (see examples section below).
##' 
##' @author Wilfried Thuiller, Damien Georges, Robin Engler
##' @seealso \code{\link{BIOMOD_FormatingData}},  
##'   \code{\link{BIOMOD_ModelingOptions}}, 
##'   \code{\link{BIOMOD_Projection}}
##' 
##' @keywords models
##' @keywords regression
##' @keywords nonlinear
##' @keywords multivariate
##' @keywords nonparametric
##' @keywords tree
##'   
##' @examples 
##' ##' species occurrences
##' DataSpecies <- 
##'   read.csv(
##'     system.file(
##'       "external/species/mammals_table.csv",
##'       package="biomod2"
##'     )
##'   )
##' head(DataSpecies)
##' 
##' ##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' ##' the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' ##' the XY coordinates of species data
##' myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##' 
##' 
##' ##' Environmental variables extracted from BIOCLIM (bio_3, 
##' ##' bio_4, bio_7, bio_11 & bio_12)
##' myExpl <- 
##'   raster::stack(
##'     system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio12.grd", package = "biomod2")
##'   )
##'
##' ##' 1. Formatting Data
##' myBiomodData <- 
##'   BIOMOD_FormatingData(
##'     resp.var = myResp,
##'     expl.var = myExpl,
##'     resp.xy = myRespXY,
##'     resp.name = myRespName
##'   )
##' 
##' ##' 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' ##' 3. Doing Modelisation
##' myBiomodModelOut <- 
##'   BIOMOD_Modeling(
##'     myBiomodData,
##'     models = c('SRE','RF'),
##'     models.options = myBiomodOption,
##'     NbRunEval = 2,
##'     DataSplit = 80,
##'     VarImport = 0,
##'     models.eval.meth = c('TSS','ROC'),
##'     do.full.models = FALSE,
##'     modeling.id = "test"
##'   )
##' 
##' ##' print a summary of modeling stuff
##' myBiomodModelOut
##' 
BIOMOD_Modeling <- function(
  data,
  models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Phillips.2'),
  models.options = NULL,
  NbRunEval = 1,
  DataSplit = 100,
  Yweights = NULL,
  Prevalence = NULL,
  VarImport = 0,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = TRUE,
  modeling.id = as.character(format(Sys.time(), "%s")),
  ...){

  # 0. loading required libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  .Models.dependencies(silent=TRUE, models.options=models.options )

  # 1. args checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  args <- .Models.check.args(data, models, models.options, NbRunEval, DataSplit,
                             Yweights, VarImport, models.eval.meth, Prevalence,
                             do.full.models, SaveObj,...)
  # updating Models arguments
  models <- args$models
  models.options <- args$models.options
  NbRunEval <- args$NbRunEval
  DataSplit <- args$DataSplit
  Yweights <- args$Yweights
  VarImport <- args$VarImport
  models.eval.meth <- args$models.eval.meth
  Prevalence <- args$Prevalence
  do.full.models <- args$do.full.models
  DataSplitTable <- args$DataSplitTable
  SaveObj <- args$SaveObj
  compress.arg = TRUE#ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))

  rm(args)
  models.out <- new('BIOMOD.models.out',
                    sp.name = data@sp.name,
                    modeling.id = modeling.id,
                    expl.var.names = colnames(data@data.env.var),
                    has.evaluation.data = data@has.data.eval,
                    rescal.all.models = rescal.all.models)

#   #To keep track of Yweights state at origin (user's input)
#     if(NbRepPA!=0 && is.null(Yweights)) Yweights <- matrix(NA, nc=Biomod.material$NbSpecies, nr=nrow(DataBIOMOD))


  # 2. creating simulation directories =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  # create the directories in which various objects will be stored (models, predictions and
  # projections). Projections' directories are created in the Projection() function.
  .Models.prepare.workdir(data@sp.name, models.out@modeling.id)


  # 3. Saving Data and Model.option objects -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(SaveObj){
    # save Input Data
    save(data, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"formated.input.data"), compress = compress.arg)
    models.out@formated.input.data@inMemory <- FALSE
    models.out@formated.input.data@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"formated.input.data")
    # save Model Options
    save(models.options, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.options"), compress = compress.arg)
    models.out@models.options@inMemory <- FALSE
    models.out@models.options@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.options")

  }


  # 3. rearanging data and determining calib and eval data -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #)
  # browser()
  mod.prep.dat <- 
    .Models.prepare.data(
      data, 
      NbRunEval, 
      DataSplit, 
      Yweights, 
      Prevalence, 
      do.full.models, 
      DataSplitTable
    )
  rm(data)

  # keeping calibLines
  calib.lines <- mod.prep.dat[[1]]$calibLines
  if(length(mod.prep.dat) > 1){ ## stack calib lines matrix along 3rd dimention of an array
    for(pa in 2:length(mod.prep.dat)){
      calib.lines <- abind(calib.lines, mod.prep.dat[[pa]]$calibLines, along=3)
    }
#     ## update dimnames
#     dimnames(calib.lines) <- list(dimnames(calib.lines)[[1]], dimnames(calib.lines)[[2]], paste("PA", 1:length(mod.prep.dat) ))
  } # else { ## force calib.line object to be a 3D array
#     dim(calib.lines) <- c(dim(calib.lines),1)
#     ## update dimnames
#     dimnames(calib.lines) <- list(dimnames(calib.lines)[[1]], dimnames(calib.lines)[[2]], paste("PA", 1:length(mod.prep.dat) ))
#   }
  # save calib.lines
  save(calib.lines, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"calib.lines"), compress = compress.arg)
  models.out@calib.lines@inMemory <- FALSE
  models.out@calib.lines@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"calib.lines")
  rm(calib.lines)


  # 4. Print modelling summary in console -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  .Models.print.modeling.summary(mod.prep.dat, models)

  # 5. Running models -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

  # loop on PA
  modeling.out <- lapply(mod.prep.dat,.Biomod.Models.loop,
                          modeling.id = models.out@modeling.id,
                          Model = models,
                          Options = models.options,
                          VarImport = VarImport,
                          mod.eval.method = models.eval.meth,
                          SavePred = SaveObj,
                          scal.models = rescal.all.models
                          )

  # put outputs in good format and save those
  # browser()
  models.out@models.computed <- .transform.outputs.list(modeling.out, out='models.run')
  models.out@models.failed <- .transform.outputs.list(modeling.out, out='calib.failure')

  if(SaveObj){
    # save model evaluation
    models.evaluation <- .transform.outputs.list(modeling.out, out='evaluation')
    save(models.evaluation, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.evaluation"), compress = compress.arg)
    models.out@models.evaluation@inMemory <- TRUE
    models.out@models.evaluation@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.evaluation")
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)

    # save model variables importances
    if(VarImport > 0 ){
      variables.importances <- .transform.outputs.list(modeling.out, out='var.import')

      ## trick to put appropriate dimnames
#       vi.dim.names <- dimnames(variables.importances)
#       vi.dim.names[[1]] <- models.out@expl.var.names
#       dimnames(variables.importances) <- vi.dim.names
#       rm('vi.dim.names')

      save(variables.importances, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"variables.importance"), compress = compress.arg)
      models.out@variables.importances@inMemory <- TRUE
      models.out@variables.importances@link <-file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"variables.importance")
      models.out@variables.importances@val <- variables.importances
      rm(variables.importances)
    }

    # save model predictions
    models.prediction <- .transform.outputs.list(modeling.out, out='prediction')
    save(models.prediction, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction"),  compress=compress.arg)
    models.out@models.prediction@inMemory <- FALSE
    models.out@models.prediction@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction")
#     models.out@models.prediction@val <- .transform.outputs(modeling.out, out='prediction')
    rm(models.prediction)

    # save evaluation model predictions
    models.prediction.eval <- .transform.outputs.list(modeling.out, out='prediction.eval')
    save(models.prediction.eval, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction.eval"), compress = compress.arg)
    models.out@models.prediction.eval@inMemory <- FALSE
    models.out@models.prediction.eval@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction.eval")
#     models.out@models.prediction@val <- .transform.outputs(modeling.out, out='prediction')
    rm(models.prediction.eval)

  }

  # removing MAXENT.Phillips tmp dir
#   if('MAXENT.Phillips' %in% models){
#     .Delete.Maxent.WorkDir(species.name=models.out@sp.name, modeling.id=models.out@modeling.id)
#   }

  rm(modeling.out)

  # save model object on hard drive
  models.out@link <- file.path(models.out@sp.name, paste(models.out@sp.name, '.', models.out@modeling.id, '.models.out', sep=""))
  assign(x=paste(models.out@sp.name, '.', models.out@modeling.id, '.models.out', sep=""),
         value=models.out)
  save(list=paste(models.out@sp.name, '.', models.out@modeling.id, '.models.out', sep=""),
       file=models.out@link)


  .bmCat("Done")
  return(models.out)
}

# -=-=-=- Several hidden functions -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #


.Models.dependencies <- function(silent=TRUE, models.options = NULL){
  # Loading all required libraries
  cat('\n\nLoading required library...')
#   require(nnet, quietly=silent)
#   require(rpart, quietly=silent)
#   require(MASS, quietly=silent)
#   require(gbm, quietly=silent)
#   require(mda, quietly=silent)
#   require(randomForest, quietly=silent)
#
#   if(!is.null(models.options)){
#     if(grepl('mgcv', models.options@GAM$algo)){
#       if("package:gam" %in% search() ) detach(package:gam)
#       require(mgcv, quietly=silent)
#     } else{
#       if("package:mgcv" %in% search() ) detach(package:mgcv)
#       require(gam, quietly=silent)
#     }
#   } else {
#     if('mgcv' %in% rownames(installed.packages())){
#       if("package:gam" %in% search() ) detach(package:gam)
#       require(mgcv, quietly=silent)
#     } else{
#       if("package:mgcv" %in% search() ) detach(package:mgcv)
#       require(gam, quietly=silent)
#     }
#   }
#
#   require(abind, quietly=silent)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.check.args <- function(data, models, models.options, NbRunEval, DataSplit,
                               Yweights, VarImport, models.eval.meth, Prevalence, do.full.models, SaveObj, ...){
  cat('\n\nChecking Models arguments...\n')
  add.args <- list(...)

  # data checking
  if(
    !inherits(
      data, 
      c("BIOMOD.formated.data", "BIOMOD.formated.data.PA", "BIOMOD.formated.data.indep", 
        "BIOMOD.formated.data.PA.indep")
    )
  ){
    stop("data argument must be a 'BIOMOD.formated.data' (obtained by running Initial.State function) ")
  }

  # models checking
  if( !is.character( models ) )
  {
    stop("models argument must be a 'character' vector")
  }

  models <- unique(models)
  models.swich.off <- NULL

  avail.models.list <- 
    c(
      'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips', 
      'MAXENT.Phillips.2'
    )
  
  ## check if model is supported
  purrr::map(models, ~ checkmate::assert_choice(.x, avail.models.list))
  # if(sum(models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Tsuruoka')) != length(models)){
  #   stop(paste(models[which( (models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Tsuruoka'))
  #                            == FALSE) ]," is not an availabe model !",sep=""))
  # }

  categorial_var <- unlist(sapply(colnames(data@data.env.var), function(x){if(is.factor(data@data.env.var[,x])) return(x) else return(NULL)} ))

  if(length(categorial_var)){
    models.fact.unsuprort <- c("SRE", "MAXENT.Tsuruoka")
    models.swich.off <- c(models.swich.off, intersect(models, models.fact.unsuprort))
    if(length(models.swich.off)){
      models <- setdiff(models, models.swich.off)
      cat(paste0("\n\t! ", paste(models.swich.off, collapse = ",", sep = " ")," were switched off because of categorical variables !"))
    }
  }
  
  ## disable MAXENT.Tsuruoka because of package maintaining issue (request from B Ripley 03-2019)
  if('MAXENT.Tsuruoka' %in% models){
    models.swich.off <- unique(c(models.swich.off, "MAXENT.Tsuruoka"))
    models <- setdiff(models, models.swich.off)
    warning('MAXENT.Tsuruoka has been disabled because of package maintaining issue (request from cran team 03-2019)')
  }

  # models.options checking ( peut etre permetre l'utilisation de liste de params )
  if(!is.null(models.options) & !inherits(models.options, "BIOMOD.Model.Options")){
    stop("models.options argument must be a 'BIOMOD.Model.Options.object' (obtained by running ... ) ")
  }

  if(is.null(models.options)){
    warning("Models will run with 'defaults' parameters", immediate.=T)
    # create a default models.options object
    models.options <- BIOMOD_ModelingOptions() # MAXENT.Phillips = list( path_to_maxent.jar = getwd())

  }

  # MAXENT.Phillips specific checking
  if("MAXENT.Phillips" %in% models){
    if(!file.exists(file.path(models.options@MAXENT.Phillips$path_to_maxent.jar ,"maxent.jar")) ){
      models = models[-which(models=='MAXENT.Phillips')]
      warning("The maxent.jar file is missing. You need to download this file (http://www.cs.princeton.edu/~schapire/maxent) and put the maxent.jar file in your working directory -> MAXENT.Phillips was switched off")
    ## -- 
    ## The java installation check is temporally disable cause it seems to cause 
    ## issues on some Windows users machine.
    ## --
    # } else if(!.check.java.installed()){
    #   models = models[-which(models=='MAXENT.Phillips')]
    } else if(nrow(data@coord)==1){
     # no coordinates
      warning("You must give XY coordinates if you want to run MAXENT.Phillips -> MAXENT.Phillips was switched off")
      models = models[-which(models=='MAXENT.Phillips')]
    }
  }

  ## Data split checks
  if(!is.null(add.args$DataSplitTable)){
    cat("\n! User defined data-split table was given -> NbRunEval, DataSplit and do.full.models argument will be ignored")
    if(!(length(dim(add.args$DataSplitTable)) %in% c(2,3) )) stop("DataSplitTable must be a matrix or a 3D array")
    if(dim(add.args$DataSplitTable)[1] != length(data@data.species) ) stop("DataSplitTable must have as many rows (dim1) than your species as data")
    NbRunEval <- dim(add.args$DataSplitTable)[2]
    DataSplit <- 50
    do.full.models <- FALSE
  }


  # NbRunEval checking (permetre un nb different par espece?)
  if( !is.numeric(NbRunEval) || NbRunEval <= 0 ){
    stop("NbRunEval argument mus be a non null positive 'numeric'")
  }

  # DataSplit checking
  if( !is.numeric(DataSplit) || DataSplit < 0 || DataSplit > 100 ){
    stop("DataSplit argument must be a 0-100 'numeric'")
  }

  if(DataSplit < 50){
    warning("You chose to allocate more data to evaluation than to calibration of your model
            (DataSplit<50)\nMake sure you really wanted to do that. \n", immediate.=T)
  }

#   # EM weight checking
#   if(!is.null(EM.weight))
#     if(!any(EM.weight==c("Roc","TSS","Kappa")))
#       stop("The 'EM.weight' parameter must be one of the following: NULL, 'Roc', 'TSS' or 'Kappa'.\n")
#
  # Check that the weight matrix was entered correctly
  if(!is.null(Yweights)){
     if(!is.numeric(Yweights))
        stop("Yweights must be a numeric vector")
     if(length(Yweights) != length(data@data.species))
       stop("The number of 'Weight' does not match with the input calibration data.
            Simulation cannot proceed.")
  }

  # Defining evaluation runs.
  if(NbRunEval <= 0){
      DataSplit <- 100
      if(!inherits(data, c("BIOMOD.formated.data.indep", "BIOMOD.formated.data.PA.indep"))){
        warning("The models will be evaluated on the calibration data only (NbRunEval=0 and no
                independent data) \n\t it could lead to over-optimistic predictive performances.\n",
                immediate.=T)
      }
  }
  if(DataSplit==100) NbRunEval <- 0

  # Models evaluation method checking
  models.eval.meth <- unique(models.eval.meth)

  if(sum(models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS','KAPPA','ACCURACY','BIAS',
                              'POD','PODFD','CSI','ETS','HK','ROC')) != length(models.eval.meth)){
    stop(paste(models.eval.meth[which( (models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS',
                                                                'KAPPA','ACCURACY','BIAS', 'POD',
                                                                'PODFD','CSI', 'ETS','HK','ROC'))
                                       == FALSE) ]," is not a availabe models evaluation metric !",sep=""))
  }

  # Prevalence checking
  if(!is.null(Prevalence)){
    if(!is.numeric(Prevalence) | Prevalence>=1 | Prevalence <=0){
      stop("Prevalence must be a 0-1 numeric")
    } else {
      # update MAXENT.Phillips default prevalence
      if("MAXENT.Phillips" %in% models){
        cat("\n\t MAXENT.Phillips defaultprevalence option was updated to fit with modeling prevalence (i.e",Prevalence,")")
        models.options@MAXENT.Phillips$defaultprevalence = Prevalence
      }
    }
  }

  ##### TO BE CHANGE BUT PREVENT FROM BUGS LATTER
  # Force object saving parameter
  if(!SaveObj){
    cat("\n\t SaveObj param was automatically set to TRUE to prevent bugs.")
    SaveObj <- TRUE
  }

#   cat('\nChecking done!\n')
  return(list(models = models,
              models.options = models.options,
              NbRunEval = NbRunEval,
              DataSplit = DataSplit,
              Yweights = Yweights,
              VarImport = VarImport,
              models.eval.meth = models.eval.meth,
              Prevalence = Prevalence,
              do.full.models = do.full.models,
              SaveObj = SaveObj,
              DataSplitTable=add.args$DataSplitTable))
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.prepare.workdir <- function(sp.name, modeling.id){
  cat("\nCreating suitable Workdir...\n")
  dir.create(sp.name, showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(sp.name,".BIOMOD_DATA",modeling.id), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(sp.name, "models",modeling.id), showWarnings=FALSE, recursive=T)

#   if(sum(models.list %in% c('MARS', 'FDA', 'ANN')) > 0 ){
#     dir.create(paste(getwd(),"/",sp.name, "/models/scaling_models", sep=""), showWarnings=FALSE, recursive=T)
#   }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.SampleMat <- function(data.sp, dataSplit, nbRun = 1, data.env = NULL){
  # return a matrix with nbRun columns of boolean (T: calib, F= eval)
  # data.sp is a 0,1 vector
  pres <- which(data.sp == 1)
  abs <- (1:length(data.sp))[-pres]

  nbPresEval <- round(length(pres) * dataSplit/100)
  nbAbsEval <- round(length(abs) * dataSplit/100)

  mat.out <- matrix(FALSE,
                    nrow = length(data.sp),
                    ncol = nbRun)
  colnames(mat.out) <- paste('_RUN',1:nbRun, sep='')

  for (i in 1:ncol(mat.out)){
    ## force to sample at least one level of each factorial variable for calibration
    fact.cell.samp <- NULL
    if(!is.null(data.env)){
      fact.cell.samp <- sample.factor.levels(data.env)
      mat.out[fact.cell.samp, i] <- TRUE
    }
    mat.out[sample(setdiff(pres, fact.cell.samp),
                   max(nbPresEval - length(fact.cell.samp), 0)), i] <- TRUE
    mat.out[sample(setdiff(abs, fact.cell.samp),
                   max(nbAbsEval - length(fact.cell.samp), 0)), i] <- TRUE
  }
  return(mat.out)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.print.modeling.summary <- function( mod.prep.dat, models){
  cat("\n\n")
  .bmCat(paste(unlist(strsplit(mod.prep.dat[[1]]$name,'_'))[1], "Modeling Summary"))

  cat("\n",ncol(mod.prep.dat[[1]]$dataBM)-1, " environmental variables (", colnames(mod.prep.dat[[1]]$dataBM)[-1], ")")

  cat("\nNumber of evaluation repetitions :" , ncol(mod.prep.dat[[1]]$calibLines))

  cat("\nModels selected :", models, "\n")

  cat("\nTotal number of model runs :",ncol(mod.prep.dat[[1]]$calibLines) * length(models) * length(mod.prep.dat),"\n")

  .bmCat()
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.check.EF.args <- function(models, models.eval.meth, models.options){
  # the models selected for doing EF
  if(models.options@EF$models.selected == 'all'){
    EF.algo <- models[which(models != 'SRE')] # remove SRE technique if it was selected (SRE cannot be used for ensemble forecast)
  } else {
    EF.algo <- models[models %in% models.options@EF$models.selected]
  }
  if(length(EF.algo)==0) stop('No models available selected for Ensemble forecasting stuff')
  # the weight methods
  if(models.options@EF$weight.method == 'all'){
    EF.weight <- models.eval.meth
  } else {
    EF.weight <- models.eval.meth[models.eval.meth %in% models.options@EF$weight.method]
  }
  if(length(EF.weight)==0) stop('No weighting method available selected for Ensemble forecasting stuff')
  return(list(EF.algo = EF.algo,
            EF.weight = EF.weight))
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.check.java.installed <- function(){

  if(.Platform$OS.type == "unix"){
    java.test <- try( expr = eval(system("command -v java", intern=TRUE )) , silent=TRUE)
  } else if(.Platform$OS.type == "windows"){
    java.test <- try( expr = eval(system( "java", intern=TRUE )), silent=TRUE )
  } else java.test <- ""

  if(!is.null(attr(java.test,"class"))){
    cat("\n! java software seems not be correctly installed\n  > MAXENT.Phillips modelling was switched off!")
    return(FALSE)
  } else{ return(TRUE) }

}

#####################################################################################################
#' Reshape biomod2 objects
#' 
#' This is an internal function (developper only)
#'
#' @param modOut the model object to transform given as a list
#' @param out character, the type of output to be transformed
#'
#' @return extracted statistics of interest from the model object
#'   as `array`.
#' @export
#'
.transform.outputs.array <-
  function(
    modOut, 
    out = 'evaluation'
  ){
    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'var.import', 'calib.failure', 'models.run', 'prediction.eval' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'var.import',
                                                            'calib.failure', 'models.run', 'prediction.eval' ))))
    }
    
    # check dim of input list
    if(length(dim(modOut)) != 4 ){
      cat('\n',dim(modOut),'\n')
      print(dimnames(modOut))
      warning("Not computed .transform.outputs because of an incompatible input list dimension", immediate=T)
      return(NULL)
    }
    
    if(dim(modOut)[4] == 1 & length(unlist(strsplit(unlist(dimnames(modOut)[4]),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(length(dimnames(modOut)[[4]]) > 0){
        dataset.names <- unlist(sapply(unlist(dimnames(modOut)[4]), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else {
        dataset.names <- paste('PA', 1:dim(modOut)[4])
      }
    }
    
    run.eval.names <- sub('_','',unlist(dimnames(modOut)[3]))
    mod.names <- unlist(dimnames(modOut)[2])
    
    if (out=='evaluation'){
      if( is.null(modOut['evaluation',1,1,1])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut['evaluation',1,1,1]))
      eval.col.names <- colnames(as.data.frame(modOut['evaluation',1,1,1]))
      
      eval.out <- array(data = unlist(modOut['evaluation',,,]),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                        eval.col.names,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(eval.out)
    }
    
    if (out=='prediction'){
      if( is.null(modOut['pred',1,1,1])){ return(NULL) }
      nb.pts.pred <- length(as.numeric(unlist(modOut['pred',1,1,1])))
      pred.out <- array(data = unlist(modOut['pred',,,]),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(pred.out)
    }
    
    if (out=='prediction.eval'){
      if( is.null(modOut['pred.eval',1,1,1])){ return(NULL) }
      nb.pts.pred.eval <- length(as.numeric(unlist(modOut['pred.eval',1,1,1])))
      pred.eval.out <- array(data = unlist(modOut['pred.eval',,,]),
                             dim = c(nb.pts.pred.eval,
                                     length(mod.names),
                                     length(run.eval.names),
                                     length(dataset.names)),
                             dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
      return(pred.eval.out)
    }
    
    if (out=='var.import'){
      if( is.null(unlist(modOut['var.import',1,1,1]))){ return(NULL) }
      nb.var <- length(as.numeric(unlist(modOut['var.import',1,1,1])))
      
      vi.out <- array(data = unlist(modOut['var.import',,,]),
                      dim = c(nb.var,
                              length(mod.names),
                              length(run.eval.names),
                              length(dataset.names)),
                      dimnames = list(paste('Var',1:nb.var,sep=''), # to change
                                      mod.names,
                                      run.eval.names,
                                      dataset.names))
      
      return(vi.out)
    }
    
    if (out == 'calib.failure'){
      cf.out <- unlist(modOut['calib.failure',,,])
      return(cf.out[!is.null(cf.out)])
    }
    
    if (out == 'models.run'){
      mod.run.out <- unlist(modOut['ModelName',,,])
      return(mod.run.out[!is.null(mod.run.out)])
    }
    
  }

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #


#' Reshape biomod2 objects
#' 
#' This is an internal function (developper only)
#'
#' @param modOut the object to transform given as a list
#' @param out character, the type of input object
#' @param dim.names character, if not `NULL` the resshaped object will be stored on the hard drive
#'
#' @return
#' list, the extracted statistics
#' @export
.transform.outputs.list =
  function(
    modOut, 
    out = 'evaluation', 
    dim.names = NULL
  ){
    
    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
                    'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'prediction.eval', 'var.import',
                                                            'calib.failure', 'models.run', 'EF.prediction',
                                                            'EF.PCA.median', 'EF.evaluation'))))
    }
    
    if(length(modOut) == 1 & length(unlist(strsplit(unlist(names(modOut)),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(is.null(dim.names)){
        dataset.names <- unlist(sapply(unlist(names(modOut)), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else{
        dataset.names <- unlist(dim.names[1])
      }
    }
    
    if(is.null(dim.names)){
      run.eval.names <- sub('_','',unlist(names(modOut[[1]]))) # may be good here to test that all names are identics
      
      mod.names <- unlist(names(modOut[[1]][[1]]))
    } else{
      run.eval.names <- unlist(dim.names[2])
      mod.names <- unlist(dim.names[3])
    }
    
    if (out=='evaluation'){
      
      eval.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            eval.tab <- modOut[[i]][[j]][[k]][['evaluation']]
            if(!is.null(eval.tab)){ break }
          }
          if(!is.null(eval.tab)){ break }
        }
        if(!is.null(eval.tab)){ break }
      }
      
      if( is.null(eval.tab)){ return(NULL) }
      
      eval.meth.names <- rownames(as.data.frame(eval.tab))
      eval.col.names <- colnames(as.data.frame(eval.tab))
      
      eval.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            if(is.null(modOut[[d1]][[d2]][[d3]][['calib.failure']])){
              return(data.frame(modOut[[d1]][[d2]][[d3]][['evaluation']]))
            } else { matrix(NA, ncol=length(eval.col.names), nrow=length(eval.meth.names), dimnames=list(eval.meth.names,eval.col.names))}
          })
        })
      })
      
      eval.out <- array(data = unlist(eval.out),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                        eval.col.names,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(eval.out)
    }
    
    if (out=='prediction'){
      
      pred.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            pred.tab <- modOut[[i]][[j]][[k]][['pred']]
            if(!is.null(pred.tab)){ break }
          }
          if(!is.null(pred.tab)){ break }
        }
        if(!is.null(pred.tab)){ break }
      }
      
      if( is.null(pred.tab)){ return(NULL) }
      
      
      nb.pts.pred <- length(as.numeric(pred.tab))
      
      pred.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            if(is.null(modOut[[d1]][[d2]][[d3]][['pred']])){
              return(rep(NA,nb.pts.pred))
            } else{
              return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred']]))
            }
          })
        })
      })
      
      pred.out <- array(data = unlist(pred.out),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(pred.out)
    }
    
    if (out=='prediction.eval'){
      pred.eval.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            pred.eval.tab <- modOut[[i]][[j]][[k]][['pred.eval']]
            if(!is.null(pred.eval.tab)){ break }
          }
          if(!is.null(pred.eval.tab)){ break }
        }
        if(!is.null(pred.eval.tab)){ break }
      }
      
      if( is.null(pred.eval.tab)){ return(NULL) }
      
      
      nb.pts.pred.eval <- length(as.numeric(pred.eval.tab))
      
      pred.eval.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            if(is.null(modOut[[d1]][[d2]][[d3]][['pred.eval']])){
              return(rep(NA,nb.pts.pred.eval))
            } else{
              return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred.eval']]))
            }
          })
        })
      })
      
      pred.eval.out <- array(data = unlist(pred.eval.out),
                             dim = c(nb.pts.pred.eval,
                                     length(mod.names),
                                     length(run.eval.names),
                                     length(dataset.names)),
                             dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
      return(pred.eval.out)
    }
    
    if (out=='var.import'){
      vi.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            vi.tab <- modOut[[i]][[j]][[k]][['var.import']]
            if(!is.null(vi.tab)){ break }
          }
          if(!is.null(vi.tab)){ break }
        }
        if(!is.null(vi.tab)){ break }
      }
      
      if( is.null(vi.tab)){ return(NULL) }
      
      nb.var <- length(as.numeric(unlist(vi.tab)))
      
      ef.mod <- grep(pattern="EF.",mod.names) # EF models
      if(length(ef.mod)>0){
        kept.mod <- mod.names[-ef.mod]
      } else{
        kept.mod <- mod.names
      }
      
      vi.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(kept.mod, function(d3){ # models without EF ones
            if(is.null(modOut[[d1]][[d2]][[d3]][['var.import']])){
              return(rep(NA,nb.var))
            } else{
              return(as.matrix(modOut[[d1]][[d2]][[d3]][['var.import']]))
            }
          })
        })
      })
      
      vi.out <- array(data = unlist(vi.out),
                      dim = c(nb.var,
                              length(kept.mod),
                              length(run.eval.names),
                              length(dataset.names)),
                      dimnames = list(names(modOut[[1]][[1]][[1]][['var.import']]), # to change
                                      kept.mod,
                                      run.eval.names,
                                      dataset.names))
      
      return(vi.out)
    }
    
    if (out == 'calib.failure'){
      cf.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            return(modOut[[d1]][[d2]][[d3]][['calib.failure']])
          })
        })
      })
      cf.out <- unlist(cf.out)
      if(length(cf.out)) cf.out <- na.omit(cf.out)
      if(length(cf.out)) cf.out <- cf.out[!is.null(cf.out)]
      if(!length(cf.out)) cf.out <- 'none'
      return(cf.out)
    }
    
    if (out == 'models.run'){
      mod.run.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            return(as.character(modOut[[d1]][[d2]][[d3]][['ModelName']]))
          })
        })
      })
      mod.run.out <- unlist(mod.run.out)
      if(length(mod.run.out)) mod.run.out <- na.omit(mod.run.out)
      if(length(mod.run.out)) mod.run.out <- mod.run.out[!is.null(mod.run.out)]
      if(!length(mod.run.out)) mod.run.out <- 'none'
      return(mod.run.out)
    }
    
    
    if (out == 'EF.prediction'){
      if( is.null(modOut[[1]][[1]][[1]][['EM']])){ return(NULL) }
      
      nb.pts.ef.pred <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['EM']])))
      
      ef.pred.out <- lapply(1:length(modOut),function(d1){ # data set
        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
            return(as.numeric(modOut[[d1]][[d2]][[d3]][['EM']]))
          })
        })
      })
      
      ef.pred.out <- array( data = unlist(ef.pred.out),
                            dim = c(nb.pts.ef.pred,
                                    length(modOut[[1]][[1]]),
                                    length(modOut[[1]]),
                                    length(modOut)),
                            dimnames = list(NULL,
                                            mod.names,
                                            run.eval.names,
                                            dataset.names))
      
      return(ef.pred.out)
    }
    
    if (out == 'EF.PCA.median'){
      if( is.null(modOut[[1]][[1]][[1]][['PCA.median']])){ return(NULL) }
      
      ef.pca.out <- lapply(1:length(modOut),function(d1){ # data set
        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
            return(as.character(modOut[[d1]][[d2]][[d3]][['PCA.median']]))
          })
        })
      })
      
      ef.pca.out <- array( data = unlist(ef.pca.out),
                           dim = c(1,
                                   length(modOut[[1]][[1]]),
                                   length(modOut[[1]]),
                                   length(modOut)),
                           dimnames = list(NULL,
                                           mod.names,
                                           run.eval.names,
                                           dataset.names))
      
      return(ef.pca.out)
    }
    
    if (out == 'EF.evaluation'){
      if( is.null(modOut[[1]][[1]][[1]][['EM.eval']])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      
      ef.eval.out <- lapply(1:length(modOut),function(d1){ # data set
        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
            return(data.frame(modOut[[d1]][[d2]][[d3]][['EM.eval']]))
          })
        })
      })
      
      ef.eval.out <- array(data = unlist(ef.eval.out),
                           dim = c(length(eval.meth.names),
                                   length(eval.col.names),
                                   length(modOut[[1]][[1]]),
                                   length(modOut[[1]]),
                                   length(modOut)),
                           dimnames = list(eval.meth.names,
                                           eval.col.names,
                                           mod.names,
                                           run.eval.names,
                                           dataset.names))
      
      return(ef.eval.out)
    }
    
  }

DF_to_ARRAY <- function(df){
  if(!is.data.frame(df) & !is.matrix(df)){
    if(is.list(df)){
      df.names <- names(df)
      df <- as.data.frame(df)
      names(df) <- df.names
    } else{
      stop("You have to give a data.frame")
    }
  }
  
  a <- sapply(strsplit(colnames(df), '_'), tail, n=3)
  b <- lapply(1:3, function(id) return(unique(a[id,])))
  array.dim.names <- c(list(character(0)),rev(b))
  #   array.dim.names <- c(list(c(NULL)),rev(apply(sapply(strsplit(colnames(df), '_'), tail, n=3),1,unique)))
  
  array.dim <- c(nrow(df),sapply(array.dim.names[-1],length))
  array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
  
  for(x in colnames(df)){
    dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
    array.out[,dimTmp[1],dimTmp[2],dimTmp[3]] <- df[,x]
  }
  return(array.out)
}

LIST_to_ARRAY <- function(ll){
  test <- sapply(ll, is.array)
  if(!all(test)) stop("list elements should be arrays")
  test <- sapply(ll,dim)
  test <- apply(test,1,function(x){length(unique(x))==1})
  if(!all(test)) stop("list elements differ in dimension")
  
  formal.dim.names <- dimnames(ll[[1]])
  new.dim.names <- rev(apply(sapply(strsplit(names(ll), '_'), tail, n=3),1,unique))
  array.dim.names <- c(formal.dim.names,new.dim.names)
  array.dim <- sapply(array.dim.names,length)
  
  array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
  
  for(x in names(ll)){
    dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
    dimTmp <- paste( paste(rep(",",length(formal.dim.names)),collapse="") , paste("'",dimTmp,"'",sep="",collapse=","),collapse="")
    
    eval(parse(text=paste("array.out[",dimTmp,"] <-  ll[[x]]",sep="")))
  }
  return(array.out)
}


##' @name BIOMOD_ModelingOptions
##' @aliases BIOMOD_ModelingOptions
##' @title Configure the modeling options for each selected model
##'
##' @description Parametrize and/or tune biomod's single models options.
##'
##' @usage
##'   BIOMOD_ModelingOptions(GLM = NULL,
##'                          GBM = NULL,
##'                          GAM = NULL,
##'                          CTA = NULL,
##'                          ANN = NULL,
##'                          SRE = NULL,
##'                          FDA = NULL,
##'                          MARS = NULL,
##'                          RF = NULL,
##'                          MAXENT.Phillips = NULL)
##'
##' @param GLM  list, GLM options
##' @param GBM  list, GBM options
##' @param GAM  list, GAM options
##' @param CTA  list, CTA options
##' @param ANN  list, ANN options
##' @param SRE  list, SRE options
##' @param FDA  list, FDA options
##' @param MARS list, MARS options
##' @param RF list, RF options
##' @param MAXENT.Phillips  list, MAXENT.Phillips options
##'
##'
##' @details
##'   The aim of this function is to allow advanced user to change some default parameters of BIOMOD inner models.
##'   For each modeling technique, options can be set up.
##'
##'   Each argument have to be put in a list object.
##'
##'   The best way to use this function is to print defaut models options (\code{\link{Print_Default_ModelingOptions}}) or create a default 'BIOMOD.model.option object' and print it in your console. Then copy the output, change only the required parameters, and paste it as function arguments. (see example)
##'
##'   Here the detailed list of modifiable parameters. They correspond to the traditional parameters that could be setted out for each modeling technique (e.g. ?GLM)
##'
##' @section GLM (\code{\link[stats]{glm}}):
##'
##'   \itemize{
##'
##'     \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off.
##'       You can choose to either:
##'         \itemize{
##'           \item{generate automatically the GLM formula by using the type and interaction.level arguments
##'             type (default \code{'quadratic'}) : formula given to the model ('simple', 'quadratic' or 'polynomial').
##'             interaction.level (default \code{0}) : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GLM.}
##'           \item{or construct specific formula}
##'         }}
##'
##'     \item{\code{test} (default \code{'AIC'}) : Information criteria for the stepwise selection procedure: AIC for Akaike Information Criteria, and BIC for Bayesian Information Criteria ('AIC' or 'BIC'). 'none' is also a supported value which implies to concider only the full model (no stepwise selection). This can lead to convergence issu and strange results.}
##'
##'     \item{\code{family} (default \code{binomial(link = 'logit')}) : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default.}
##'
##'     \item{\code{control} : a list of parameters for controlling the fitting process. For glm.fit this is passed to \code{\link{glm.control}}.}
##'
##'   }
##'
##' @section GBM (default \code{\link[gbm]{gbm}}):
##'
##'   Please refer to \code{\link[gbm]{gbm}} help file to get the meaning of this options.
##'   \itemize{
##'     \item{ \code{distribution} (default \code{'bernoulli'})}
##'     \item{ \code{n.trees} (default \code{2500})}
##'     \item{ \code{interaction.depth} (default \code{7})}
##'     \item{ \code{n.minobsinnode} (default \code{5})}
##'     \item{ \code{shrinkage} (default \code{0.001})}
##'     \item{ \code{bag.fraction} (default \code{0.5})}
##'     \item{ \code{train.fraction} (default \code{1})}
##'     \item{ \code{cv.folds} (default \code{3})}
##'     \item{ \code{keep.data} (default \code{FALSE})}
##'     \item{ \code{verbose} (default \code{FALSE})}
##'     \item{ \code{perf.method} (default \code{'cv'})}
##'     \item{ \code{n.cores} (default \code{1})}
##'   }
##'
##'
##' @section GAM (\code{\link[gam]{gam}} or \code{\link[mgcv]{gam}}):
##'   \itemize{
##'
##'     \item{algo : either "GAM_gam" (default), "GAM_mgcv" or "BAM_mgcv" defining the chosen GAM function (see \code{\link[mgcv]{gam}}, \code{\link[gam]{gam}} resp. \code{\link[mgcv]{bam}} for more details)}
##'
##'     \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off.
##'       You can choose to either:
##'         \itemize{
##'           \item{generate automatically the GAM formula by using the type and interaction.level arguments
##'             type : the smother used to generate the formula. Only "s_smoother" available at time.
##'             interaction.level : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GAM. Interaction are not considered if you choosed "GAM_gam" algo}
##'           \item{or construct specific formula}
##'         }}
##'
##'     \item{k (default \code{-1} or \code{4}): a smooth term in a formula argument to gam (see \pkg{gam} \code{\link[gam]{s}} or \pkg{mgcv} \code{\link[mgcv]{s}})}
##'     \item{family (default \code{binomial(link = 'logit')}) : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default. }
##'     \item{control : see \code{\link[mgcv]{gam.control}} or \code{\link[gam]{gam.control}}}
##'
##'     \item{some extra "GAM_mgcv" specific options (ignored if algo = "GAM_gam")
##'       \itemize{
##'         \item{\code{method} (default \code{'GCV.Cp'})}
##'         \item{\code{optimizer} (default \code{c('outer','newton')})}
##'         \item{\code{select} (default \code{FALSE})}
##'         \item{\code{knots} (default \code{NULL})}
##'         \item{\code{paramPen} (default \code{NULL})}
##'       }
##'
##'     }
##'   }
##'
##'
##' @section CTA (\code{\link[rpart]{rpart}}):
##'
##'   Please refer to \code{\link[rpart]{rpart}} help file to get the meaning of the following options.
##'   \itemize{
##'     \item{\code{method} (default \code{'class'})}
##'     \item{\code{parms} (default \code{'default'}) : if \code{'default'}, default \pkg{rpart} parms value are kept}
##'     \item{\code{cost} (default \code{NULL})}
##'     \item{\code{control}: see \code{\link[rpart]{rpart.control}}}
##'   }
##'
##'   NOTE: for method and parms, you can give a 'real' value as described in the rpart help file or 'default' that implies default \code{\link[rpart]{rpart}} values.
##'
##' @section ANN (\code{\link[nnet]{nnet}}):
##'
##'   \itemize{
##'     \item{\code{NbCV} (default \code{5}) : nb of cross validation to find best size and decay parameters}
##'     \item{\code{size}} (default \code{NULL}) : number of units in the hidden layer. If \code{NULL} then size parameter will be optimised by cross validation based on model AUC (\code{NbCv} cross validation; tested size will be the following c(2,4,6, 8) ). You can also specified a vector of size you want to test. The one giving the best model AUC will be then selected.
##'     \item{\code{decay}} (default \code{NULL}) : parameter for weight decay. If \code{NULL} then decay parameter will be optimised by cross validation on model AUC (\code{NbCv} cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1) ). You can also specified a vector of decay you want to test. The one giving the best model AUC will be then selected.
##'     \item{\code{rang} (default \code{0.1}) : Initial random weights on [-rang, rang]}
##'     \item{\code{maxit} (default \code{200}): maximum number of iterations.}
##'   }
##'
##'
##' @section SRE (\code{\link[biomod2]{sre}}):
##'   \itemize{
##'     \item{\code{quant} (default \code{0.025}): quantile of 'extreme environmental variable' removed for selection of species envelops}
##'   }
##'
##'
##' @section FDA (\code{\link[mda]{fda}}):
##'
##'   Please refer to \code{\link[mda]{fda}} help file to get the meaning of these options.
##'   \itemize{
##'     \item{\code{method} (default \code{'mars'})}
##'     \item{\code{add_args} (default \code{NULL}) : additional arguments to \code{method} given as a list of parameters
##'       (corespond to the \ldots options of fda function) }
##'   }
##'
##' @section MARS (\code{\link[earth]{earth}}):
##'
##'   Please refer to \code{\link[earth]{earth}} help file to get the meaning of these options.
##'   \itemize{
##'     \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off.
##'       You can choose to either:
##'         \itemize{
##'           \item{generate automatically the GLM formula by using the type and interaction.level arguments
##'             type (default \code{'simple'}) : formula given to the model ('simple', 'quadratic' or 'polynomial').
##'             interaction.level (default \code{0}) : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GLM/MARS.}
##'           \item{or construct specific formula}
##'         }}
##'     %    \item{\code{degree} (default \code{2})}
##'     \item{\code{nk}} (default \code{NULL}) : an optional integer specifying the maximum number of model terms. If NULL is given then default mars function value is used ( i.e max(21, 2 * nb_expl_var + 1) )
##'     \item{\code{penalty} (default \code{2})}
##'     \item{\code{thresh} (default \code{0.001})}
##'     \item{\code{nprune} (default \code{NULL})}
##'     \item{\code{pmethod} (default \code{"backward"})}
##'   }
##'
##'
##' @section RF (\code{\link[randomForest]{randomForest}}):
##'
##'   \itemize{
##'     \item{\code{do.classif} (default \code{TRUE}) : if TRUE classification random.forest computed else regression random.forest will be done}
##'     \item{\code{ntree} (default \code{500})}
##'     \item{\code{mtry} (default \code{'default'})}
##'     \item{\code{nodesize} (default \code{5})}
##'     \item{\code{maxnodes} (default \code{NULL})}
##'   }
##'
##'   NOTE: for mtry, you can give a 'real' value as described in randomForest help file or 'default' that implies default randomForest values
##'
##' @section  [MAXENT.Phillips](https://biodiversityinformatics.amnh.org/open_source/maxent/) :
##'   \itemize{
##'     \item{\code{path_to_maxent.jar} : character, the link to \pkg{maxent.jar} file (the working directory by default) }
##'     \item{\code{memory_allocated} : integer (default \code{512}), the amount of memory (in Mo) reserved for java to run MAXENT.Phillips. should be 64, 128, 256, 512, 1024, 2048... or NULL if you want to use default java memory limitation parameter.}
##'     \item{\code{background_data_dir} : character, path to a directory where explanatory variables are stored as ASCII files (raster format).
##'       If specified MAXENT will generate it's own background data from expalantory variables rasters (as usually done in MAXENT studies). If not
##'       set, then MAXENT will use the same pseudo absences than other models (generated within biomod2 at formatting step) as background data.}
##'     \item{\code{maximumbackground} : integer, the maximum number of background data to sample. This parameter will be use only if \code{background_data_dir}
##'       option has been set to a non default value.}
##'     \item{\code{maximumiterations} : integer (default \code{200}), maximum iteration done}
##'     \item{\code{visible} : logical (default \code{FALSE}), make the Maxent user interface visible}
##'     \item{\code{linear} : logical (default \code{TRUE}), allow linear features to be used}
##'     \item{\code{quadratic} : logical (default \code{TRUE}), allow quadratic features to be used}
##'     \item{\code{product} : logical (default \code{TRUE}), allow product features to be used}
##'     \item{\code{threshold} : logical (default \code{TRUE}), allow threshold features to be used}
##'     \item{\code{hinge} : logical (default \code{TRUE}), allow hinge features to be used}
##'     \item{\code{lq2lqptthreshold} : integer (default \code{80}), number of samples at which product and threshold features start being used}
##'     \item{\code{l2lqthreshold} : integer (default \code{10}), number of samples at which quadratic features start being used}
##'     \item{\code{hingethreshold} : integer (default \code{15}), number of samples at which hinge features start being used}
##'     \item{\code{beta_threshold} : numeric (default \code{-1.0}), regularization parameter to be applied to all threshold features; negative value enables automatic setting}
##'     \item{\code{beta_categorical} : numeric (default \code{-1.0}), regularization parameter to be applied to all categorical features; negative value enables automatic setting}
##'     \item{\code{beta_lqp} : numeric (default \code{-1.0}), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting}
##'     \item{\code{beta_hinge} : numeric (default \code{-1.0}), regularization parameter to be applied to all hinge features; negative value enables automatic setting}
##'     \item{\code{betamultiplier} : numeric (default \code{1}), multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution.}
##'     \item{\code{defaultprevalence} : numeric (default \code{0.5}), default prevalence of the species: probability of presence at ordinary occurrence points}
##'   }
##'
##' % @section MAXENT.Tsuruoka (\code{\link[maxent]{maxent}}):
##' %
##' % \itemize{
##' %   \item{\code{l1_regularizer} (default \code{0.0}): An numeric turning on L1 regularization and setting the regularization parameter. A value of 0 will disable L1 regularization}
##' %   \item{\code{l2_regularizer} (default \code{0.0}): An numeric turning on L2 regularization and setting the regularization parameter. A value of 0 will disable L2 regularization}
##' %   \item{\code{use_sgd} (default \code{FALSE}): A logical indicating that SGD parameter estimation should be used. Defaults to FALSE}
##' %   \item{\code{set_heldout} (default \code{0}): An integer specifying the number of documents to hold out. Sets a held-out subset of your data to test against and prevent overfitting}
##' %   \item{\code{verbose} (default \code{FALSE}): A logical specifying whether to provide descriptive output about the training process}
##' % }
##'
##' % NOTE: if you use the \code{set_heldout} parameter then the data that will be held out will be taken in the
##' % calibration data pool. It can be penilizing in case of low number of occurences dataset.
##'
##'
##' @return
##'   A \code{"\link[=BIOMOD.Model.Options-class]{BIOMOD.Model.Options}"} object given to \code{\link[biomod2]{BIOMOD_Modeling}}
##'
##' @author Damien Georges, Wilfried Thuiller
##' @keywords models
##' @keywords options
##'
##' @examples
##'   ## default BIOMOD.model.option object
##'   myBiomodOptions <- BIOMOD_ModelingOptions()
##'
##'   ## print the object
##'   myBiomodOptions
##'
##'   ## you can copy a part of the print, change it and custom your options
##'   ## here we want to compute quadratic GLM and select best model with 'BIC' criterium
##'   myBiomodOptions <- BIOMOD_ModelingOptions(
##'     GLM = list( type = 'quadratic',
##'                 interaction.level = 0,
##'                 myFormula = NULL,
##'                 test = 'BIC',
##'                 family = 'binomial',
##'                 control = glm.control(epsilon = 1e-08,
##'                                       maxit = 1000,
##'                                       trace = FALSE) ))
##'
##'   ## check changes was done
##'   myBiomodOptions
##'
##'   ##' you can prefer to establish your own GLM formula
##'   myBiomodOptions <- BIOMOD_ModelingOptions(
##'     GLM = list( myFormula = formula("Sp277 ~ bio3 +
##'                     log(bio10) + poly(bio16,2) + bio19 + bio3:bio19")))
##'
##'   ## check changes was done
##'   myBiomodOptions
##'
##'   ##' you also can directly print default parameters and then follow the same processus
##'   Print_Default_ModelingOptions()
##'


####################################################################################################
'BIOMOD_ModelingOptions' <- function(
                        GLM = NULL,
                        GBM = NULL,
                        GAM = NULL,
                        CTA = NULL,
                        ANN = NULL,
                        SRE = NULL,
                        FDA = NULL,
                        MARS = NULL,
                        RF = NULL,
                        MAXENT.Phillips = NULL
                        ){
  # 1. create a defaut BIOMOD.Model.Options object
  opt <- new('BIOMOD.Model.Options')

  # 2. modify it if necessary
  if(!is.null(GLM)){
    if(!is.null(GLM$type)) { opt@GLM$type <- GLM$type }
    if(!is.null(GLM$interaction.level)) { opt@GLM$interaction.level <- GLM$interaction.level }
    if(!is.null(GLM$myFormula)) { opt@GLM$myFormula <- GLM$myFormula }
    if(!is.null(GLM$test)) { opt@GLM$test <- GLM$test }
    if(!is.null(GLM$family)) {
      fam.test <- TRUE
      if(inherits(GLM$family, 'family')){
        opt@GLM$family <- GLM$family
      } else{
        if( is.character(GLM$family)){
          if(! unlist(strsplit(GLM$family,"[/(]"))[1] %in% c('binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson')){ fam.test <- FALSE}

          if(grepl(')', GLM$family)){ # check string formalisation to add () if necessary
            opt@GLM$family <- eval(parse(text=GLM$family))
          } else{
            opt@GLM$family <- eval(parse(text=paste(GLM$family,"()", sep="")))
          }
        } else{ fam.test <- FALSE }
      }
      if(!fam.test){
        cat("\n!!! invalid GLM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GLM$family <- binomial(link = 'logit')
      }
    }
    if(!is.null(GLM$mustart)) { opt@GLM$mustart <- GLM$mustart }
    if(!is.null(GLM$control)) { opt@GLM$control <- GLM$control }
  }

  if(!is.null(GBM)){
#     if(!is.null(GBM$type )) { opt@GBM$type <- GBM$type }
#     if(!is.null(GBM$interaction.level )) { opt@GBM$interaction.level <- GBM$interaction.level }
    if(!is.null(GBM$distribution )) { opt@GBM$distribution <- GBM$distribution }
    if(!is.null(GBM$n.trees )) { opt@GBM$n.trees <- GBM$n.trees }
    if(!is.null(GBM$interaction.depth )) { opt@GBM$interaction.depth <- GBM$interaction.depth }
    if(!is.null(GBM$n.minobsinnode )) { opt@GBM$n.minobsinnode <- GBM$n.minobsinnode }
    if(!is.null(GBM$shrinkage )) { opt@GBM$shrinkage <- GBM$shrinkage }
    if(!is.null(GBM$bag.fraction )) { opt@GBM$bag.fraction <- GBM$bag.fraction }
    if(!is.null(GBM$train.fraction )) { opt@GBM$train.fraction <- GBM$train.fraction }
    if(!is.null(GBM$cv.folds )) { opt@GBM$cv.folds <- GBM$cv.folds }
    if(!is.null(GBM$keep.data )) { opt@GBM$keep.data <- GBM$keep.data }
    if(!is.null(GBM$verbose )) { opt@GBM$verbose <- GBM$verbose }
#     if(!is.null(GBM$class.stratify.cv )) { opt@GBM$class.stratify.cv <- GBM$cv.folds }
    if(!is.null(GBM$perf.method )) { opt@GBM$perf.method <- GBM$perf.method }
    if(!is.null(GBM$n.cores)) { opt@GBM$n.cores <- GBM$n.cores } else { opt@GBM$n.cores <- NULL }
  }



  if(!is.null(GAM)){
    if(!is.null(GAM$algo )) { opt@GAM$algo <- GAM$algo }
    if(!is.null(GAM$type )) { opt@GAM$type <- GAM$type }
    if(!is.null(GAM$k )) { opt@GAM$k <- GAM$k } else{
      if(opt@GAM$algo == 'GAM_gam'){
        opt@GAM$k <- 4
      } else{
        opt@GAM$k <- -1
      }
    }
    if(!is.null(GAM$interaction.level )) { opt@GAM$interaction.level <- GAM$interaction.level }
    if(!is.null(GAM$myFormula )) { opt@GAM$myFormula <- GAM$myFormula }
    if(!is.null(GAM$family)) {
      fam.test <- TRUE
      if(inherits(GAM$family, 'family')){
        opt@GAM$family <- GAM$family
      } else{
        if( is.character(GAM$family)){
          if(! unlist(strsplit(GAM$family,"[/(]"))[1] %in% c('binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson')){ fam.test <- FALSE}

          if(grepl(')', GAM$family)){ # check string formalisation to add () if necessary
            opt@GAM$family <- eval(parse(text=GAM$family))
          } else{
            opt@GAM$family <- eval(parse(text=paste(GAM$family,"()", sep="")))
          }
        } else{ fam.test <- FALSE }
      }
      if(!fam.test){
        cat("\n!!! invalid GAM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GAM$family <- binomial(link = 'logit')
      }
    }

    if(is.null(GAM$control )) {
      if(opt@GAM$algo == 'GAM_gam'){
        requireNamespace('gam', quietly = TRUE)
        opt@GAM$control <- gam::gam.control()
      } else{ opt@GAM$control <- mgcv::gam.control() }
    } else{
      user.control.list <- GAM$control
      if(opt@GAM$algo == 'GAM_gam'){
        default.control.list <- gam::gam.control()
      } else{
        default.control.list <- mgcv::gam.control()
      }

      control.list <- lapply(names(default.control.list), function(x){
        if(x %in% names(user.control.list)){
          return(user.control.list[[x]])
        } else {
          return(default.control.list[[x]])
        }
      })

      names(control.list) <- names(default.control.list)
      opt@GAM$control <- control.list
    }

    if(!is.null(GAM$method )) { opt@GAM$method <- GAM$method }
    if(!is.null(GAM$optimizer )) { opt@GAM$optimizer <- GAM$optimizer }
    if(!is.null(GAM$select )) { opt@GAM$select <- GAM$select }
    if(!is.null(GAM$knots )) { opt@GAM$knots <- GAM$knots }
    if(!is.null(GAM$paraPen )) { opt@GAM$paraPen <- GAM$paraPen }
  } else{
    if(opt@GAM$algo == 'GAM_gam'){
      opt@GAM$control <- gam::gam.control()
      opt@GAM$k <- 4
    } else{
      opt@GAM$control <- mgcv::gam.control()
      opt@GAM$k <- -1
    }
  }


  if(!is.null(CTA)){
#     if(!is.null(CTA$type )) { opt@CTA$type <- CTA$type }
#     if(!is.null(CTA$interaction.level )) { opt@CTA$interaction.level <- CTA$interaction.level }
    if(!is.null(CTA$method )) { opt@CTA$method <- CTA$method }
    if(!is.null(CTA$parms )) { opt@CTA$parms <- CTA$parms }
    if(!is.null(CTA$control )) { opt@CTA$control <- CTA$control }
    if(!is.null(CTA$cost )) { opt@CTA$cost <- CTA$cost }
  }

  if(!is.null(ANN)){
#     if(!is.null(ANN$type )) { opt@ANN$type <- ANN$type }
#     if(!is.null(ANN$interaction.level )) { opt@ANN$interaction.level <- ANN$interaction.level }
    if(!is.null(ANN$NbCV )) { opt@ANN$NbCV <- ANN$NbCV }
    if(!is.null(ANN$size )) { opt@ANN$size <- ANN$size }
    if(!is.null(ANN$decay )) { opt@ANN$decay <- ANN$decay }
    if(!is.null(ANN$rang )) { opt@ANN$rang <- ANN$rang }
    if(!is.null(ANN$maxit )) { opt@ANN$maxit <- ANN$maxit }
  }

  if(!is.null(SRE)){
    if(!is.null(SRE$quant )) { opt@SRE$quant <- SRE$quant }
  }

  if(!is.null(FDA)){
#     if(!is.null(FDA$type )) { opt@FDA$type <- FDA$type }
#     if(!is.null(FDA$interaction.level )) { opt@FDA$interaction.level <- FDA$interaction.level }
    if(!is.null(FDA$method )) { opt@FDA$method <- FDA$method }
    if(!is.null(FDA$add_args )) { opt@FDA$add_args <- FDA$add_args } ## additional args such as degree, nk
  }

  if(!is.null(MARS)){
    if(!is.null(MARS$type)) { opt@MARS$type <- MARS$type }
    if(!is.null(MARS$interaction.level)) { opt@MARS$interaction.level <- MARS$interaction.level }
    if(!is.null(MARS$myFormula)) { opt@MARS$myFormula <- MARS$myFormula }
#     if(!is.null(MARS$degree )) { opt@MARS$degree <- MARS$degree }
    if(!is.null(MARS$nk )) { opt@MARS$nk <- MARS$nk }
    if(!is.null(MARS$penalty )) { opt@MARS$penalty <- MARS$penalty }
    if(!is.null(MARS$thresh )) { opt@MARS$thresh <- MARS$thresh }
    if(!is.null(MARS$nprune )) { opt@MARS$nprune <- MARS$nprune }
    if(!is.null(MARS$pmethod )) { opt@MARS$pmethod <- MARS$pmethod }
  }

  if(!is.null(RF)){
    if(!is.null(RF$type )) { opt@RF$type <- RF$type }
#     if(!is.null(RF$interaction.level )) { opt@RF$interaction.level <- RF$interaction.level }
#     if(!is.null(RF$do.classif )) { opt@RF$do.classif <- RF$do.classif }
    if(!is.null(RF$ntree )) { opt@RF$ntree <- RF$ntree }
    if(!is.null(RF$mtry )) { opt@RF$mtry <- RF$mtry }
    if(!is.null(RF$nodesize )) { opt@RF$nodesize <- RF$nodesize }
    if(!is.null(RF$maxnodes )) { opt@RF$maxnodes <- RF$maxnodes }
  }

  if(!is.null(MAXENT.Phillips)){
    if(!is.null(MAXENT.Phillips$path_to_maxent.jar )) {
      opt@MAXENT.Phillips$path_to_maxent.jar <- normalizePath(sub("maxent.jar", "", MAXENT.Phillips$path_to_maxent.jar)) # ensure path format validity
      } else {opt@MAXENT.Phillips$path_to_maxent.jar <- getwd()}
    if(!is.null(MAXENT.Phillips$memory_allocated )) { opt@MAXENT.Phillips$memory_allocated <- MAXENT.Phillips$memory_allocated }
	if(!is.null(MAXENT.Phillips$background_data_dir )) { opt@MAXENT.Phillips$background_data_dir <- MAXENT.Phillips$background_data_dir }
    if(!is.null(MAXENT.Phillips$maximumbackground )) { opt@MAXENT.Phillips$maximumbackground <- MAXENT.Phillips$maximumbackground }
    if(!is.null(MAXENT.Phillips$maximumiterations )) { opt@MAXENT.Phillips$maximumiterations <- MAXENT.Phillips$maximumiterations }
    if(!is.null(MAXENT.Phillips$visible )) { opt@MAXENT.Phillips$visible <- MAXENT.Phillips$visible }
    if(!is.null(MAXENT.Phillips$linear )) { opt@MAXENT.Phillips$linear <- MAXENT.Phillips$linear }
    if(!is.null(MAXENT.Phillips$quadratic )) { opt@MAXENT.Phillips$quadratic <- MAXENT.Phillips$quadratic }
    if(!is.null(MAXENT.Phillips$product )) { opt@MAXENT.Phillips$product <- MAXENT.Phillips$product }
    if(!is.null(MAXENT.Phillips$threshold )) { opt@MAXENT.Phillips$threshold <- MAXENT.Phillips$threshold }
    if(!is.null(MAXENT.Phillips$hinge )) { opt@MAXENT.Phillips$hinge <- MAXENT.Phillips$hinge }
    if(!is.null(MAXENT.Phillips$lq2lqptthreshold )) { opt@MAXENT.Phillips$lq2lqptthreshold <- MAXENT.Phillips$lq2lqptthreshold }
    if(!is.null(MAXENT.Phillips$l2lqthreshold )) { opt@MAXENT.Phillips$l2lqthreshold <- MAXENT.Phillips$l2lqthreshold }
    if(!is.null(MAXENT.Phillips$hingethreshold )) { opt@MAXENT.Phillips$hingethreshold <- MAXENT.Phillips$hingethreshold }
    if(!is.null(MAXENT.Phillips$beta_threshold )) { opt@MAXENT.Phillips$beta_threshold <- MAXENT.Phillips$beta_threshold }
    if(!is.null(MAXENT.Phillips$beta_categorical )) { opt@MAXENT.Phillips$beta_categorical <- MAXENT.Phillips$beta_categorical }
    if(!is.null(MAXENT.Phillips$beta_lqp )) { opt@MAXENT.Phillips$beta_lqp <- MAXENT.Phillips$beta_lqp }
    if(!is.null(MAXENT.Phillips$beta_hinge )) { opt@MAXENT.Phillips$beta_hinge <- MAXENT.Phillips$beta_hinge }
	  if(!is.null(MAXENT.Phillips$betamultiplier )) { opt@MAXENT.Phillips$betamultiplier <- MAXENT.Phillips$betamultiplier }
    if(!is.null(MAXENT.Phillips$defaultprevalence )) { opt@MAXENT.Phillips$defaultprevalence <- MAXENT.Phillips$defaultprevalence }
  } else{
    opt@MAXENT.Phillips$path_to_maxent.jar <- getwd()
  }

  # if(!is.null(MAXENT.Tsuruoka)){
  #   if(!is.null(MAXENT.Tsuruoka$l1_regularizer )) { opt@MAXENT.Tsuruoka$l1_regularizer <- MAXENT.Tsuruoka$l1_regularizer }
  #   if(!is.null(MAXENT.Tsuruoka$l2_regularizer )) { opt@MAXENT.Tsuruoka$l2_regularizer <- MAXENT.Tsuruoka$l2_regularizer }
  #   if(!is.null(MAXENT.Tsuruoka$use_sgd )) { opt@MAXENT.Tsuruoka$use_sgd <- MAXENT.Tsuruoka$use_sgd }
  #   if(!is.null(MAXENT.Tsuruoka$set_heldout )) { opt@MAXENT.Tsuruoka$set_heldout <- MAXENT.Tsuruoka$set_heldout }
  #   if(!is.null(MAXENT.Tsuruoka$verbose )) { opt@MAXENT.Tsuruoka$verbose <- MAXENT.Tsuruoka$verbose }
  # }

  test <- as.logical(validObject(object = opt, test = TRUE, complete = FALSE))

  if(!test){
    cat("\n\n!!! NULL object returned because of invalid parameters given !!!")
    return(NULL)
  }

  return(opt)
}

Print_Default_ModelingOptions <- function(){
  cat('\n Defaut modeling options. copy, change what you want paste it as arg to BIOMOD_ModelingOptions\n\n')

  opt_tmp <- BIOMOD_ModelingOptions()
  print(opt_tmp)
}

BIOMOD_presenceonly <- function(modeling.output = NULL, EM.output = NULL, bg.env = NULL, perc = 0.9, save.output = T){
  # if(!require(PresenceAbsence)){stop("PresenceAbsence package required!")}
  requireNamespace('PresenceAbsence', quietly = TRUE)
  
  myModelEval <- myBiomodProjFF <- NULL
  
  if(!is.null(modeling.output)){  
    calib.lines<-get(load(modeling.output@calib.lines@link))[,,1]
    myResp <- get(load(modeling.output@formated.input.data@link))@data.species
    
    myModelEval <- get_evaluations(modeling.output,as.data.frame=T)
    myModelEval[,1] <- as.character(myModelEval[,1])
    for(i in 1:nrow(myModelEval)){myModelEval[i,1] <- paste(c(modeling.output@sp.name,strsplit(as.character(myModelEval[i,1]),split="_")[[1]][3:1]),collapse="_")  } 
    
    myModelPred <- get_predictions(modeling.output,as.data.frame=T)
    if(!is.null(bg.env)){
      myModelPred.pres <- myModelPred[myResp==1,]
      myBiomodProj.eval <- BIOMOD_Projection(
        new.env = bg.env,
        proj.name = paste(modeling.output@modeling.id,"cv_EF_eval",sep="_"),
        modeling.output = modeling.output,
        build.clamping.mask = F)      
      myModelPred <- as.data.frame(myBiomodProj.eval@proj@val)
      ### Change the colnames to the real model names
      colnames(myModelPred) <- 
        paste(
          modeling.output@sp.name,
          rep(dimnames(myBiomodProj.eval@proj@val)[[4]],prod(dim(myBiomodProj.eval@proj@val)[2:3])),
          rep(dimnames(myBiomodProj.eval@proj@val)[[3]],each = dim(myBiomodProj.eval@proj@val)[2]),
          rep(dimnames(myBiomodProj.eval@proj@val)[[2]],dim(myBiomodProj.eval@proj@val)[3]), sep='_')
    }
    if(modeling.output@has.evaluation.data == T){
      myModelPred.eval  <- as.data.frame(get(load(paste(modeling.output@"sp.name","/.BIOMOD_DATA/",modeling.output@modeling.id,"/models.prediction.eval", sep=""))))
      for(i in 1:ncol(myModelPred.eval)){colnames(myModelPred.eval)[i] <- paste(c(modeling.output@sp.name,strsplit(colnames(myModelPred.eval)[i],split="[.]")[[1]][3:1]),collapse="_")  }       
    }
  }
  
  if(!is.null(EM.output)){
    if(EM.output@em.by!='PA_dataset+repet'){stop("em.by of 'BIOMOD.EnsembleModeling' must be 'PA_dataset+repet'")} 
    myModelEvalEF <- get_evaluations(EM.output,as.data.frame=T)
    myModelEvalEF[,1] <- paste(modeling.output@sp.name,as.character(myModelEvalEF[,1]),sep="_")
    #for(i in 1:nrow(myModelEvalEF)){myModelEvalEF[i,1] <- paste("EF",strsplit(as.character(myModelEvalEF[,1]),split="_")[[i]][2],"AllData",sep="_")}
    if(!is.null(modeling.output)){
      myModelEval <- rbind(myModelEval, myModelEvalEF)
    }
    
    myBiomodProjFF <- get_predictions(EM.output,as.data.frame=T)  
    
    if(!is.null(bg.env)){
      myBiomodProjFF.pres <- as.data.frame(myBiomodProjFF[myResp==1,])
      colnames(myBiomodProjFF.pres) <- colnames(myBiomodProjFF)
      myBiomodProjFF <- BIOMOD_EnsembleForecasting(
        proj.name = paste(modeling.output@modeling.id,"cv_EF_bg",sep="_"), 
        projection.output = myBiomodProj.eval,
        EM.output = EM.output)    
      myBiomodProjFF <- as.data.frame(myBiomodProjFF@proj@val)     
      myModelPred.pres <- cbind(myModelPred.pres,myBiomodProjFF.pres)
    }
    myModelPred <- cbind(myModelPred, myBiomodProjFF)
    
    if(modeling.output@has.evaluation.data == T){
      myBiomodProjFF.eval <- get_predictions(EM.output,as.data.frame=T,evaluation=T)  
      
      #colnames(myBiomodProjFF.eval) <- gsub("AllAlgos_ROC_EMwmean","EF",  colnames(myBiomodProjFF.eval))
      myModelPred.eval <- cbind(myModelPred.eval, myBiomodProjFF.eval)      
    }  
  }
  
  mpa.eval <- boyce.eval <- myModelEval[!duplicated(myModelEval[,1]),]
  boyce.eval$Eval.metric <- "boyce"
  mpa.eval$Eval.metric <- "mpa"
  boyce.eval[,3:7]<-mpa.eval[,3:7]<-NA
  
  ###MPA & BOYCE     
  for(i in 1:nrow(boyce.eval)){
    n <- length(strsplit(as.character(boyce.eval[i,1]),split="_")[[1]])
    tec <- paste(strsplit(as.character(boyce.eval[i,1]),split="_")[[1]][3:n],collapse="_") 
    Model.name <- boyce.eval[i,1]
    run <- strsplit(Model.name,split="_")[[1]][c(grep("RUN",strsplit(Model.name,split="_")[[1]]),grep("Full",strsplit(Model.name,split="_")[[1]]))]
    
    #### CORRECTION ------------------------------------------------------------------------
    if(inherits(calib.lines, "matrix")){
      ind.eval = which(calib.lines[,paste("_",run, sep="")] == FALSE) #### CORRECTION
    }else{
      ind.eval = which(calib.lines == FALSE) #### CORRECTION
    }
    
    if(length(ind.eval)==0){      #this is the full model ##### CORRECTION
      if(is.null(bg.env)){
        # if(inherits(calib.lines, "matrix")){ #### PROBLEM : this part gives problem with the cbind after (not same number of rows)
        #   test <- myResp
        #   Pred<-myModelPred[,Model.name]
        # }else{
          test <- myResp[calib.lines]     
          Pred <- myModelPred[calib.lines,Model.name]
        # }
      }else{
        test <- c(myResp[myResp==1],rep(0,nrow(bg.env)))    
        Pred <- c(myModelPred.pres[,Model.name],myModelPred[,Model.name])       
      }
    }else{
      if(is.null(bg.env)){
        test <- myResp[ind.eval] #### CORRECTION
        Pred <- myModelPred[ind.eval,Model.name] #### CORRECTION
      }else{
        test <- c(myResp[ind.eval & myResp==1],rep(0,nrow(bg.env))) #### CORRECTION
        Pred <- c(myModelPred.pres[ind.eval & myResp==1,Model.name],myModelPred[,Model.name])  #### CORRECTION
      }
    }
    #### CORRECTION ------------------------------------------------------------------------
    
    
    ind.1 = which(test == 1)
    ind.notNA = which(!is.na(Pred))
    ind.obs = intersect(ind.1, ind.notNA)
    
    if (length(ind.obs) > 0){ #### CORRECTION
      boy <- ecospat::ecospat.boyce(fit=Pred[ind.notNA],obs=Pred[ind.obs], PEplot=F) #### CORRECTION
      boyce.eval[boyce.eval[,1]==Model.name,3] <- boy$Spearman.cor
      if( sum(boy$F.ratio<1,na.rm=T)>0){
        boyce.eval[boyce.eval[,1]==Model.name,5] <- round(boy$HS[max(which(boy$F.ratio<1))],0)
        DATA<-cbind(1:length(Pred), test, Pred/1000) #### PROBLEM
        DATA[is.na(DATA[,2]),2] <- 0
        DATA <- DATA[stats::complete.cases(DATA),]
        if(!is.na(round(boy$HS[max(which(boy$F.ratio<1))],0)/1000)){
          EVAL<-presence.absence.accuracy(DATA, threshold=round(boy$HS[max(which(boy$F.ratio<1))],0)/1000) 
          boyce.eval[boyce.eval[,1]==Model.name,6] <-  EVAL$sensitivity 
          boyce.eval[boyce.eval[,1]==Model.name,7] <-  EVAL$specificity
        }else{boyce.eval[boyce.eval[,1]==Model.name,6:7] <-  NA}
      }else{
        boyce.eval[boyce.eval[,1]==Model.name,7] <-  boyce.eval[boyce.eval[,1]==Model.name,6] <-  boyce.eval[boyce.eval[,1]==Model.name,5] <- NA 	
      }
      
      mpa.eval[mpa.eval[,1]==Model.name,5] <- ecospat::ecospat.mpa(Pred[ind.obs], perc = perc) #### CORRECTION
      EVAL<-presence.absence.accuracy(DATA, threshold=ecospat::ecospat.mpa(Pred[ind.obs], perc = perc)/1000) #### CORRECTION
      mpa.eval[mpa.eval[,1]==Model.name,6] <-  EVAL$sensitivity 
      mpa.eval[mpa.eval[,1]==Model.name,7] <-  EVAL$specificity  
    }
    
    if(modeling.output@has.evaluation.data == T){
      myResp.eval <- get(load(modeling.output@formated.input.data@link))@eval.data.species
      Pred.eval <- myModelPred.eval[,Model.name]
      
      boy <- ecospat::ecospat.boyce(fit=Pred.eval,obs=Pred.eval[myResp.eval==1 & ind.1], PEplot=F) #### CORRECTION
      boyce.eval[boyce.eval[,1]==Model.name,"Evaluating.data"] <- boy$Spearman.cor
      
      mpa.eval[mpa.eval[,1]==Model.name,"Evaluating.data"] <- ecospat::ecospat.mpa(Pred.eval[myResp.eval==1 & ind.1], perc = perc) #### CORRECTION
    }
  }
  myModelEval[,6:7] <- round(myModelEval[,6:7],1)
  boyce.eval[,6:7] <- round(boyce.eval[,6:7]*100,1)
  mpa.eval[,6:7] <- round(mpa.eval[,6:7]*100,1)  
  
  if(modeling.output@has.evaluation.data == T){
    if(!is.null(EM.output)){
      output <- list(eval=rbind(myModelEval,boyce.eval,mpa.eval),myBiomodProjFF=myBiomodProjFF,myBiomodProjEF.eval=myBiomodProjFF.eval) 
    }else{
      output <- list(eval=rbind(myModelEval,boyce.eval,mpa.eval))      
    }
    if(save.output){
      if(!is.null(modeling.output)){
        sp<-modeling.output@sp.name}
      if(!is.null(EM.output)){
        sp<-EM.output@sp.name}      
      save(output, paste(sp, "/.BIOMOD_DATA/",modeling.output@modeling.id,"/presenceonlyevaluation_",sp, sep=""))
    }
    return(output)
  }else{
    if(!is.null(EM.output)){
      output <- list(eval=rbind(myModelEval,boyce.eval,mpa.eval),myBiomodProjFF=myBiomodProjFF)
    }else{
      output <- list(eval=rbind(myModelEval,boyce.eval,mpa.eval))      
    }
    if(save.output){
      if(!is.null(modeling.output)){
        sp<-modeling.output@sp.name}
      if(!is.null(EM.output)){
        sp<-EM.output@sp.name}      
      save(output, file=paste(sp, "/.BIOMOD_DATA/",modeling.output@modeling.id,"/presenceonlyevaluation_",sp, sep=""))
    }
    return(output)
  }
}
####################################################################################################
# BIOMOD_Projection
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   Project models from BIOMOD_Modeling with different explanatory variables

# INPUT :


# OUTPUT :



# NOTE :
#   It would be nice to add done projections to input Biomod.models.object
#   .BIOMOD_Projection.check.args <- may be reorder variables if necessary

####################################################################################################
'BIOMOD_Projection' <- function(modeling.output,
                                new.env,
                                proj.name,
                                xy.new.env = NULL,
                                selected.models = 'all',
                                binary.meth = NULL,
                                filtered.meth = NULL,
                                compress = TRUE,
                                build.clamping.mask = TRUE,
                                ...){

  # 0. get additional args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  args <- list(...)
  omit.na <- args$omit.na # omit all non full filled environemental cells (always TRUE if env is a raster object)
  silent <- args$silent # echo advancement or not
  do.stack <- args$do.stack # store output in a lone stack
  clamping.level <- args$clamping.levels # remove all cells where at least clamping.level variables are out of their calibrating range
#   clamping.value <- args$clamping.value # reference value for clamped cells
  output.format <- args$output.format # raster output format
  keep.in.memory <- args$keep.in.memory # store results on memory or only on hard drive
  on_0_1000 <- args$on_0_1000 # transform projections on a 0 - 1000 scale to limit memory consumption
  split.proj <- args$split.proj

  if(is.null(omit.na)) omit.na <- TRUE
  if(is.null(silent)) silent <- FALSE
  if(is.null(do.stack)) do.stack <- TRUE
#   if(is.null(clamping.value)) clamping.value <- -1
  if(is.null(keep.in.memory)) keep.in.memory <- TRUE
  if(is.null(on_0_1000)) on_0_1000 <- TRUE # by default we return projections on a 0 -  1000 scale.
  if(is.null(split.proj)) split.proj <- 1 # by default we return projections on a 0 -  1000 scale.

#   if(!do.stack | !keep.in.memory) rasterOptions(todisk=TRUE)

#   if(is.null(output.format)){
#     if(!inherits(new.env,"Raster"))
#       output.format <- ".RData"
#     else
#       output.format <- ".grd"
#   }

  if(!silent){
    .bmCat("Do Models Projections")
  }

  # 1. Some Inputs args checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  args <- .BIOMOD_Projection.check.args(modeling.output,
                                        new.env,
                                        proj.name,
                                        xy.new.env,
                                        selected.models,
                                        binary.meth,
                                        filtered.meth,
                                        compress,
                                        do.stack,
                                        output.format)#, clamping.level)

  proj.name <- args$proj.name
  selected.models <- args$selected.models
  binary.meth <- args$binary.meth
  filtered.meth <- args$filtered.meth
  compress <- args$compress
  do.stack <- args$do.stack
  xy.new.env <- args$xy.new.env
  output.format <- args$output.format
#   clamping.level <- args$clamping.level

  rm(args)

  # 1b. Creating the outpput object
  proj_out <- new('BIOMOD.projection.out',
              proj.names = proj.name,
              sp.name =  modeling.output@sp.name,
              expl.var.names = modeling.output@expl.var.names,
              models.projected = selected.models,
              scaled.models = modeling.output@rescal.all.models,
              xy.coord = xy.new.env,
              modeling.object.id = modeling.output@modeling.id)

  # add link to biomod2 modeling object
  proj_out@modeling.object@link = modeling.output@link


  # adapting the proj slot to projection data type (e.g. rasterStack, array)
#   if(!do.stack){
#     proj_out@proj <- new('BIOMOD.stored.files')
#   } else if(inherits(new.env, 'Raster')){
#     proj_out@proj <- new('BIOMOD.stored.raster.stack')
#   } else{
#     proj_out@proj <- new('BIOMOD.stored.array')
#   }
  if(inherits(new.env, 'Raster')){
    proj_out@proj <- new('BIOMOD.stored.raster.stack')
  } else{
    proj_out@proj <- new('BIOMOD.stored.array')
  }

  # 1.c creating output directory
  dir.create(file.path(modeling.output@sp.name,paste("proj_", proj.name, sep="")),
             showWarnings = FALSE, recursive = TRUE, mode = "777")

  # 1.c Define the clamping mask
  if(build.clamping.mask){
    if(!silent) cat("\n\t> Building clamping mask\n")
    MinMax <- get_formal_data(modeling.output,'MinMax')

    assign(x = paste("proj_",proj.name,"_",modeling.output@sp.name,"_ClampingMask",sep=""),
           value = .build.clamping.mask(new.env, MinMax) )

    if(output.format == '.RData'){
      save(list=paste("proj_",proj.name,"_",modeling.output@sp.name,"_ClampingMask",sep=""),
           file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_ClampingMask", output.format ,sep="") ), compress=compress )
    } else {
      writeRaster(x=get(paste("proj_",proj.name,"_",modeling.output@sp.name,"_ClampingMask",sep="")),
                  filename=file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_ClampingMask", output.format ,sep="")),
                  datatype = "INT2S", NAflag=-9999,
                  overwrite=TRUE)
    }


  }

  # 2. Making projections -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(!do.stack){
    dir.create(file.path(modeling.output@sp.name,paste("proj_", proj.name, sep=""), "individual_projections"),
               showWarnings = FALSE, recursive = TRUE, mode = "777")
  }

  saved.files <- c()

  proj <- lapply(selected.models, function(mod.name){
    cat("\n\t> Projecting",mod.name,"...")
    ## load biomod model
    BIOMOD_LoadModels(modeling.output, full.name=mod.name, as="mod")

    filename <- NULL
    if(!do.stack){
      filename <- file.path(modeling.output@sp.name, paste("proj_", proj.name, sep=""), "individual_projections", paste("proj_", proj.name, "_", mod.name,ifelse(output.format==".RData",".grd",output.format), sep="") )
    }

    pred.tmp <- predict(mod, new.env, on_0_1000=on_0_1000, filename=filename, omit.na=omit.na, split.proj = split.proj)

    if(do.stack){ # return the prediction only if stack has to be build
      return(pred.tmp)
    } else {
      return(filename)
    }

  })

  # 2b. Puting outputs in the right format =-=-=-=-=-=-=-=-=-=-=-= #
  if(!do.stack){
    saved.files <- unlist(proj)
  } else {
    if(inherits(new.env, "Raster")){
      proj <- stack(proj)
      names(proj) <- selected.models
    } else {
      proj <- as.data.frame(proj)
      names(proj) <- selected.models
      proj <- DF_to_ARRAY(proj)
    }
    if( keep.in.memory ){
      proj_out@proj@val <- proj
      proj_out@proj@inMemory <- TRUE
    }
  }

  ## save projections
  assign(x = paste("proj_",proj.name, "_", modeling.output@sp.name, sep=""),
         value = proj)

  if(output.format == '.RData'){
    saved.files <- file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name, output.format ,sep=""))
    save(list = paste("proj_",proj.name, "_", modeling.output@sp.name, sep=""),
         file = saved.files, compress=compress)
  } else {
    if(do.stack){
      saved.files <- file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name, output.format ,sep=""))
      writeRaster(x=get(paste("proj_",proj.name, "_", modeling.output@sp.name, sep="")),
                  filename=saved.files, overwrite=TRUE,
                  datatype = "INT2S", NAflag=-9999
                  # datatype = "INT2U"
                  )

    }
  }

  proj_out@type <- class(proj_out@proj@val)
  proj_out@proj@link <- saved.files


  # 3. Compute binary and filtered transformation =-=-=-=-=-=-=-=- #
  if(!is.null(binary.meth) | !is.null(filtered.meth)){
    cat("\n")
    eval.meth <- unique(c(binary.meth,filtered.meth))

    ## get all treshold
    if(inherits(new.env, "Raster")){
      thresholds <- matrix(0,nrow=length(eval.meth), ncol=length(selected.models),
                           dimnames=list(eval.meth,selected.models))
      for(mod in selected.models){
        PA.run   <- .extractModelNamesInfo(model.names=mod, info='data.set')
        eval.run <- .extractModelNamesInfo(model.names=mod, info='run.eval')
        algo.run <- .extractModelNamesInfo(model.names=mod, info='models')
        thresholds[eval.meth,mod] <- get_evaluations(modeling.output)[eval.meth,"Cutoff",algo.run,eval.run,PA.run]
        if(! on_0_1000) thresholds[eval.meth,mod]  <- thresholds[eval.meth,mod] / 1000
      }
    } else{
      thresholds <- array(0,dim=c(length(eval.meth),dim(proj)[-1]) ,
                         dimnames=c(list(eval.meth), dimnames(proj)[-1]) )
      for(mod in selected.models){
        PA.run   <- .extractModelNamesInfo(model.names=mod, info='data.set')
        eval.run <- .extractModelNamesInfo(model.names=mod, info='run.eval')
        algo.run <- .extractModelNamesInfo(model.names=mod, info='models')
        thresholds[eval.meth,algo.run,eval.run,PA.run] <- get_evaluations(modeling.output)[eval.meth,"Cutoff",algo.run,eval.run,PA.run]
        if(! on_0_1000) thresholds[eval.meth,algo.run,eval.run,PA.run]  <- thresholds[eval.meth,algo.run,eval.run,PA.run] / 1000
      }
    }

    ## do binary transformation
    for(eval.meth in binary.meth){
      cat("\n\t> Building", eval.meth,"binaries")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)[,i]
          writeRaster(x = BinaryTransformation(raster(file.tmp, RAT=FALSE),thres.tmp),
                      filename = sub(output.format, paste("_",eval.meth,"bin", output.format, sep=""), file.tmp),
                      overwrite=TRUE,
                      datatype = "INT2S",NAflag=-9999)
        }
      } else {
      assign(x = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep=""),
             value = BinaryTransformation(proj,asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)))

      if(output.format == '.RData'){
        save(list = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep=""),
             file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"bin", output.format ,sep="")), compress=compress)
      } else {
        writeRaster(x=get(paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep="")),
                    filename=file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"bin", output.format ,sep="")),
                    overwrite=TRUE,
                    datatype = "INT2S", NAflag=-9999)
      }

      rm(list=paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep=""))
      }
    }

    ## do filtered transformation
    for(eval.meth in filtered.meth){
      cat("\n\t> Building", eval.meth,"filtered")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)[,i]
          writeRaster(x = FilteringTransformation(raster(file.tmp, RAT=FALSE),thres.tmp),
                      filename = sub(output.format, paste("_",eval.meth,"filt", output.format, sep=""), file.tmp),
                      overwrite=TRUE, datatype=ifelse(on_0_1000,"INT2S","FLT4S"), NAflag=-9999
                      )
        }
      } else {
        assign(x = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep=""),
               value = FilteringTransformation(proj,asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)))

        if(output.format == '.RData'){
          save(list = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep=""),
               file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"filt", output.format ,sep="")), compress=compress)
        } else {
          writeRaster(x=get(paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep="")),
                      filename=file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"filt", output.format ,sep="")),
                      overwrite=TRUE , datatype=ifelse(on_0_1000,"INT2S","FLT4S"), NAflag=-9999
                      )
        }

        rm(list=paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep=""))
      }
    }

  }
  # End compute binary and filtered transformation -=-=-=-=-=-=-=- #

#   # 2.b a posteriori clamping...
#   #### TO DO : make an a priori claming
#   if(!is.null(clamping.level)){
#     cat("\n   > clamping projections...")
#     if(inherits(proj_out@proj@val, 'Raster')){
#       proj_out@proj@val <- raster::stack(reclassify(proj_out@proj@val * reclassify( (- 1 * (clampMask)), c(-0.5,0.5,1)), c(-Inf,0,clamped.value) ))
#     } else{
#       proj_out@proj@val[which(clampMask > 0),,,] <- clamped.value
#     }
#
#   }

  # save a copy of output object without value to be lighter
  assign(paste(modeling.output@sp.name,".", proj.name, ".projection.out", sep=""), free(proj_out))
  save(list = paste(modeling.output@sp.name,".", proj.name, ".projection.out", sep=""),
       file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep=""), paste(modeling.output@sp.name,".", proj.name, ".projection.out", sep="")))

  if(!silent) .bmCat("Done")

  # 4. Returning output
  return(proj_out)
}

####################################################################################################
### Utilities Fuctions #############################################################################
####################################################################################################

.BIOMOD_Projection.check.args <- function(modeling.output, new.env, proj.name, xy.new.env,
                                          selected.models, binary.meth, filtered.meth,
                                          compress, do.stack, output.format){#, clamping.level){
  ## modeling.output
  if(!inherits(modeling.output, 'BIOMOD.models.out')){
    stop("'modeling.output' must be the result of BIOMOD_Modeling() computation")
  }

  ## new.env
  # NOTE : may be reorder variables if necessary
  if(!inherits(new.env, c('matrix', 'data.frame', 'RasterStack'))){
    stop("'new.env' must be a matrix, a data.frame or a RasterStack")
  }
  if(inherits(new.env, 'RasterStack')){
    if(sum(!(names(new.env) %in% modeling.output@expl.var.names)) > 0 ){
      stop("'new.env' layer names don't match with explanatory variables used for buiding models")
    }
  } else{
    if(sum(!(colnames(new.env) %in% modeling.output@expl.var.names)) > 0 ){
      stop("'new.env' colnames don't match with explanatory variables used for buiding models")
    }
  }

  ## proj.name
  # The projection Name
  if(is.null(proj.name)){
    stop("\nYou must define a name for Projection Outputs")
  } else{
    dir.create(paste(modeling.output@sp.name,'/proj_',proj.name,'/',sep=''),
               showWarnings=FALSE)
  }

  ## xy.new.env
  if(!is.null(xy.new.env)  & !inherits(new.env,'Raster')){
    xy.new.env = data.matrix(xy.new.env)
    if(ncol(xy.new.env) != 2 | nrow(xy.new.env) != nrow(new.env)) stop("invalid xy coordinates argument given -- dimensions mismatch !")
  } else {
    xy.new.env = matrix()
  }

  ## selected.models
  if(selected.models[1] == 'all'){
    selected.models <- modeling.output@models.computed
  } else{
    selected.models <- intersect(selected.models, modeling.output@models.computed)
  }
  if(length(selected.models) < 1){
    stop('No models selected')
  }

  # check that given models exits
  files.check <- paste(modeling.output@sp.name,'/models/',modeling.output@modeling.id,"/",selected.models,sep='')
  not.checked.files <- c(grep('MAXENT.Phillips', files.check), grep('SRE', files.check))
  if(length(not.checked.files) > 0){files.check <- files.check[-not.checked.files]}
  missing.files <- files.check[!file.exists(files.check)]
  if( length(missing.files) > 0 ){
    stop(paste("Projection files missing : ", toString(missing.files), sep=''))
    if(length(missing.files) == length(files.check)){
      stop("Impossible to find any models, might be a problem of working directory")
    }
  }

  # The binaries  and filtering transformations
  if(!is.null(binary.meth) | !is.null(filtered.meth)){
    models.evaluation <- get_evaluations(modeling.output)
    if(is.null(models.evaluation)){
      warning("Binary and/or Filtered transformations of projection not ran because of models
              evaluation information missing")
    } else{
      available.evaluation <- unique(unlist(dimnames(models.evaluation)[1]))
      if(!is.null(binary.meth)){
        if(sum(!(binary.meth %in% available.evaluation)) > 0){
          warning(paste(toString(binary.meth[!(binary.meth %in% available.evaluation)]),
                        " Binary Transformation were switched off because no corresponding",
                        " evaluation method found ", sep=""))
          binary.meth <- binary.meth[binary.meth %in% available.evaluation]
        }
      }

      if(!is.null(filtered.meth)){
        if(sum(!(filtered.meth %in% available.evaluation)) > 0){
          warning(paste(toString(filtered.meth[!(filtered.meth %in% available.evaluation)]),
                        " Filtered Transformation were switched off because no corresponding",
                        " evaluation method found ", sep=""))
          filtered.meth <- filtered.meth[filtered.meth %in% available.evaluation]
        }
      }
    }
  }

  ## compress
  if(compress == 'xz'){
    compress <- ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }

  ## do.stack
  if(!inherits(new.env, 'RasterStack')){
    if(!do.stack) cat("\n\t\t! 'do.stack' arg is always set as TRUE for data.frame/matrix dataset")
    do.stack <- TRUE
  } else{
    if(do.stack){
      # test if there is memory enough to work with RasterStack
      test = canProcessInMemory( raster::subset(new.env,1), 2*length(selected.models) + nlayers(new.env) )
      if (!test) rasterOptions(todisk=T)
#       if (!do.stack){
#         cat("\n   ! Results will be saved as individual RasterLayers because of a lack of memory !")
#       }
    }
  }

  ## output.format ##
  if(!is.null(output.format)){
    if(! output.format %in% c(".img",".grd",".RData")){
      stop("output.format argument should be one of '.img','.grd' or '.RData'\n Note : '.img','.grd' are only available if you give environmental condition as a rasterStack object")
    }
    if( output.format %in% c(".img",".grd") & !inherits(new.env,"Raster") ){
      warning("output.format was automatically set to '.RData' because environmental conditions are not given as a raster object")
    }
  }
  ## set default values
  if(is.null(output.format)){
    if(!inherits(new.env,"Raster"))
      output.format <- ".RData"
    else
      output.format <- ".grd"
  }


#   ## clamping checking
#   if(!is.null(clamping.level)){
#     if(!is.numeric(clamping.level)){
#       stop("clamping.level must be NULL or numeric")
#     }
#
#     # limit clamping level
#     if( clamping.level > length(modeling.output@expl.var.names)){
#       cat("\n   ! clamping.level was down to", length(modeling.output@expl.var.names))
#       clamping.level <- length(modeling.output@expl.var.names)
#     }
#
#     if( clamping.level < 1){
#       cat("\n   ! clamping was swich off")
#       clamping.level <- NULL
#     }
#   }


  return(list(#modeling.output = modeling.output,
              #new.env = new.env,
              proj.name = proj.name,
              xy.new.env = xy.new.env,
              selected.models = selected.models,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth,
              compress = compress,
              do.stack = do.stack,
              output.format = output.format))#, clamping.level = clamping.level))

}

####################################################################################################
####################################################################################################

.build.clamping.mask <- function(env, MinMax){
  if(inherits(env,'Raster')){ # raster case
    env <- raster::stack(env)
    env <- raster::subset(env,names(MinMax))

    # create an empty mask
#     clamp.mask <- reclassify( raster::subset(env,1, drop=TRUE), c(-Inf,Inf,0) )
    clamp.mask <- ref.mask <- raster::subset(env,1, drop=TRUE)
    clamp.mask[!is.na(clamp.mask[])] <- 0
    ref.mask[!is.na(clamp.mask[])] <- 1

    for(e.v in names(MinMax)){

      if(!is.null(MinMax[[e.v]]$min)){ # numeric variable
        clamp.mask <- clamp.mask + ( BinaryTransformation(raster::subset(env, e.v, drop=TRUE), MinMax[[e.v]]$max ) +  (ref.mask - BinaryTransformation(raster::subset(env, e.v, drop=TRUE), MinMax[[e.v]]$min )) )

      } else if(!is.null(MinMax[[e.v]]$levels)){ # factorial variable
        clamp.mask <- clamp.mask + (raster::subset(env, e.v, drop=TRUE) %in% MinMax[[e.v]]$levels)
      }
    }

    ## fix projection system

  } else if(is.data.frame(env) | is.matrix(env) | is.numeric(env)){ # matrix and data.frame case
    env <- as.data.frame(env)

    # create an empty mask
    clamp.mask <- rep(0,nrow(env))

    for(e.v in names(MinMax)){
      if(!is.null(MinMax[[e.v]]$min)){ # numeric variable
        clamp.mask <- clamp.mask + ( BinaryTransformation(env[,e.v], MinMax[[e.v]]$max ) +
          (1 - BinaryTransformation(env[,e.v], MinMax[[e.v]]$min )) )
      } else if(!is.null(MinMax[[e.v]]$levels)){ # factorial variable
        clamp.mask <- clamp.mask + (env[,e.v] %in% MinMax[[e.v]]$levels)
      }
    }

  } else{
    stop("Unsupported env arg")
  }

  return(clamp.mask)

}
setGeneric( "BIOMOD_RangeSize",
            def = function(CurrentPred, FutureProj, ...){
              standardGeneric( "BIOMOD_RangeSize" )
            } )

# The data.frame input function......................................................................
setMethod('BIOMOD_RangeSize', signature(CurrentPred='data.frame', FutureProj='data.frame' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
            args <- .BIOMOD_RangeSize.check.args( CurrentPred, FutureProj,  SpChange.Save )


            #Function that estimate the number of pixel gain, loss, stable (present and absent) by species
            Value <- c(-2, 0, -1, 1)
            CompteurSp <- function(Data, Value)
            {
              if(is.data.frame(Data)) {
                N <- dim(Data)[2]
                Compt <- as.data.frame(matrix(0, ncol=4, nrow=dim(Data)[2]))
                i <- 1
                while(i <= N) {
                  Compt[i, 1] <- length(Data[Data[,i] == Value[1], i])
                  Compt[i, 2] <- length(Data[Data[,i] == Value[2], i])
                  Compt[i, 3] <- length(Data[Data[,i] == Value[3], i])
                  Compt[i, 4] <- length(Data[Data[,i] == Value[4], i])
                  i <- i + 1
                }
              }
              return(Compt)
            }
            Diff.By.Pixel <- as.data.frame(FutureProj - 2 * CurrentPred)
            Compt.By.Models <- as.data.frame(CompteurSp(Diff.By.Pixel, Value))
            Compt.By.Models[, 5] <- (100 * Compt.By.Models[, 1])/(Compt.By.Models[, 1] + Compt.By.Models[,3])
            Compt.By.Models[, 6] <- (100 * Compt.By.Models[, 4])/(Compt.By.Models[, 1] + Compt.By.Models[,3])
            Compt.By.Models[, 7] <- Compt.By.Models[, 6] - Compt.By.Models[, 5]
            Compt.By.Models[, 8] <- Compt.By.Models[, 1] + Compt.By.Models[,3]
            Compt.By.Models[, 9] <- Compt.By.Models[,3]
            Compt.By.Models[, 10] <- Compt.By.Models[, 4] + Compt.By.Models[,3]
            dimnames(Compt.By.Models) <- list(colnames(CurrentPred), c("Loss","Stable0", "Stable1", "Gain", "PercLoss", "PercGain",
                                                                       "SpeciesRangeChange", "CurrentRangeSize",
                                                                       "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
            Output <- c()
            Output <- list(Compt.By.Models=Compt.By.Models, Diff.By.Pixel=Diff.By.Pixel)
            invisible(Output)
          })

# The array input function......................................................................
setMethod('BIOMOD_RangeSize', signature(CurrentPred='array', FutureProj='array' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            # transform arrays into data.frame
            CurrentPred <- as.data.frame(CurrentPred)
            names(CurrentPred) <- unlist(lapply(strsplit(names(CurrentPred),".", fixed=TRUE),
                                                function(x){
                                                  return(paste( x[3], x[2], x[1],sep="_"))
                                                }))

            names(CurrentPred) <- unlist(lapply(strsplit(names(CurrentPred),".", fixed=TRUE),
                                             function(x){
                                               x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                               data.set.id <- x.rev[1]
                                               cross.valid.id <- x.rev[2]
                                               algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                               model.id <- paste(data.set.id,
                                                                 cross.valid.id,
                                                                 algo.id, sep="_")
                                               return(model.id)
                                             }))

            FutureProj <- as.data.frame(FutureProj)
            names(FutureProj) <- unlist(lapply(strsplit(names(FutureProj),".", fixed=TRUE),
                                                function(x){
                                                  x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                  data.set.id <- x.rev[1]
                                                  cross.valid.id <- x.rev[2]
                                                  algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                  model.id <- paste(data.set.id,
                                                                    cross.valid.id,
                                                                    algo.id, sep="_")
                                                  return(model.id)
                                                }))
            return(BIOMOD_RangeSize(CurrentPred, FutureProj,  SpChange.Save))
          })

# The Raster input function......................................................................

setMethod('BIOMOD_RangeSize', signature(CurrentPred='RasterStack', FutureProj='RasterStack' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
            args <- .BIOMOD_RangeSize.check.args( CurrentPred, FutureProj,  SpChange.Save )

            CBS <- matrix(ncol=10, nrow=length(CurrentPred@layers), dimnames=list(names(CurrentPred),
                                                                                  c("Loss","Stable0", "Stable1", "Gain", "PercLoss", "PercGain", "SpeciesRangeChange", "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")))


            sp.stack <- stack()
            for(i in 1:length(CurrentPred@layers)){
              #DiffByPixel
              Cur <- CurrentPred@layers[[i]]
              Fut <- FutureProj@layers[[i]]
              Ras <- Fut - (Cur + Cur)
              sp.stack <- addLayer(sp.stack, Ras)

              #ComptBySpecies
              CBS[i, 1] <- length(which(Ras[]==-2))
              CBS[i, 2] <- length(which(Ras[]== 0))
              CBS[i, 3] <- length(which(Ras[]==-1))
              CBS[i, 4] <- length(which(Ras[]== 1))

              CBS[i, 5] <- round(CBS[i,1] / (CBS[i,1]+CBS[i,3]) *100, digits=3)
              CBS[i, 6] <- round(CBS[i,4] / (CBS[i,1]+CBS[i,3]) *100, digits=3)
              CBS[i, 7] <- round((CBS[i,3]+CBS[i,4]) / (CBS[i,1]+CBS[i,3]) *100 -100, digits=3)

              CBS[i, 8] <- CBS[i,1]+CBS[i,3]
              CBS[i, 9] <- CBS[i,3]
              CBS[i, 10] <- CBS[i,3]+CBS[i,4]
            }

            if(is.null(SpChange.Save)) SpChange.Save <- "NoName"
            #     assign(paste(SpChange.Save, "_Compt.By.Species", sep=""), CBS, pos=1)
            names(sp.stack) <- rownames(CBS)
            #     assign(paste(SpChange.Save, "_Diff.By.Pixel", sep=""), sp.stack, pos=1)
            return(list(Compt.By.Models = CBS,
                        Diff.By.Pixel = sp.stack))
          })

####
setMethod('BIOMOD_RangeSize', signature(CurrentPred='RasterLayer', FutureProj='RasterStack' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){

            # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
            CBS <- matrix(ncol=10, nrow=length(FutureProj@layers), dimnames=list(names(FutureProj),
                                                                                 c("Loss","Stable0", "Stable1", "Gain", "PercLoss", "PercGain", "SpeciesRangeChange", "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")))
            sp.stack <- stack()
            for(i in 1:length(FutureProj@layers)){
              #DiffByPixel
              Cur <- CurrentPred
              Fut <- FutureProj@layers[[i]]
              Ras <- Fut - (Cur + Cur)
              sp.stack <- addLayer(sp.stack, Ras)

              #ComptBySpecies
              CBS[i, 1] <- length(which(Ras[]==-2))
              CBS[i, 2] <- length(which(Ras[]== 0))
              CBS[i, 3] <- length(which(Ras[]==-1))
              CBS[i, 4] <- length(which(Ras[]== 1))

              CBS[i, 5] <- round(CBS[i,1] / (CBS[i,1]+CBS[i,3]) *100, digits=3)
              CBS[i, 6] <- round(CBS[i,4] / (CBS[i,1]+CBS[i,3]) *100, digits=3)
              CBS[i, 7] <- round((CBS[i,3]+CBS[i,4]) / (CBS[i,1]+CBS[i,3]) *100 -100, digits=3)

              CBS[i, 8] <- CBS[i,1]+CBS[i,3]
              CBS[i, 9] <- CBS[i,3]
              CBS[i, 10] <- CBS[i,3]+CBS[i,4]
            }

            if(is.null(SpChange.Save)) SpChange.Save <- "NoName"
            #     assign(paste(SpChange.Save, "_Compt.By.Species", sep=""), CBS, pos=1)
            names(sp.stack) <- rownames(CBS)
            #     assign(paste(SpChange.Save, "_Diff.By.Pixel", sep=""), sp.stack, pos=1)
            return(list(Compt.By.Models = CBS,
                        Diff.By.Pixel = sp.stack))
          })


setMethod('BIOMOD_RangeSize', signature(CurrentPred='RasterLayer', FutureProj='RasterLayer' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            BIOMOD_RangeSize(CurrentPred = CurrentPred, FutureProj = stack(FutureProj),  SpChange.Save = SpChange.Save)
          })



####################################################################
.BIOMOD_RangeSize.check.args <- function( CurrentPred, FutureProj,  SpChange.Save ){
  # dimensions checking
  if(sum(!(dim(CurrentPred) == dim(FutureProj)) > 0)){
    stop("CurrentPred & FutureProj dimensions mismatched!")
  }

  # binary checking


}
##' @name BIOMOD_tuning
##' @aliases BIOMOD_tuning
##' 
##' @title Tune models parameters
##' @description Function to tune biomod single models parameters
##'
##' @param data            BIOMOD.formated.data object returned by BIOMOD_FormatingData
##' @param models          vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips' 
##' @param models.options  BIOMOD.models.options object returned by BIOMOD_ModelingOptions. Default: BIOMOD_ModelingOptions()
##' @param trControl       global control parameters for runing (default trainControl(method="cv",summaryFunction=twoClassSummary,classProbs=T),returnData = FALSE). for details see trainControl
##' @param ctrl.GAM        specify control parameters only for GAM (default trControl)
##' @param ctrl.GBM        specify control parameters only for GBM (default trControl)
##' @param ctrl.GLM        specify control parameters only for GLM (default trControl)
##' @param ctrl.CTA        specify control parameters only for CTA (default trControl)
##' @param ctrl.RF         specify control parameters only for RF (default trControl)
##' @param ctrl.ANN        specify control parameters only for ANN (default trControl)
##' @param ctrl.MARS       specify control parameters only for MARS (default trControl)
##' @param ctrl.FDA        specify control parameters only for FDA (default trControl)
##' @param metric          metric to select the optimal model (Default ROC). TSS (maximizing Sensitivity and Specificity) is also possible. see ?train
##' @param metric.ME       metric to select the optimal model for MAXENT.Phillips (Default: ROC). One out of Mean.AUC (or ROC), Mean.AUC.DIFF, Mean.ORmin, Mean.OR10 and AICc. see ?ENMevaluate and Muscarella et al. 2014
##' @param tuneLength      see ?train (default 30)
##' @param method.RF       which classification or regression model to use for randomForest (default: "rf"). see http://topepo.github.io/caret/Random_Forest.html
##' @param method.ANN      which classification or regression model to use for artificial neural networks (default: "avNNet"). see http://topepo.github.io/caret/Neural_Network.html
##' @param method.MARS     which classification or regression model to use for mars (default: "earth"). see http://topepo.github.io/caret/Multivariate_Adaptive_Regression_Splines.html
##' @param method.GAM      which classification or regression model to use for GAM (default: "gam"). see http://topepo.github.io/caret/Generalized_Additive_Model.html
##' @param method.GLM      which classification or regression model to use for GLM: (default: 'glmStepAIC'). see http://topepo.github.io/caret/Generalized_Linear_Model.html
##' @param type.GLM        vector of modeling types choosen among 'simple', 'quadratic', 'polynomial' or 's_smoother' (default c('simple','quadratic','polynomial','s_smoother'))
##' @param interaction.GLM vector of interaction type choosen among 0, 1. Default c(0,1)
##' @param cvmethod.ME     method used for data partitioning for MAXENT.Phillips (default: 'randomkfold')
##' @param kfolds.ME       number of bins to use for k-fold cross-validation used for MAXENT.Phillips (Default: 10).
##' @param overlap.ME      logical; Calculates pairwise metric of niche overlap if TRUE (Default: FALSE). (see ?calc.niche.overlap)
##' @param clamp.ME        logical; If TRUE (Default) "Features are constrained to remain within the range of values in the training data" (Elith et al. 2011)
##' @param n.bg.ME         Number of Background points used to run MAXENT.Phillips (Default: 10000)
##' @param env.ME          RasterStack of model predictor variables
##' @param size.tune.ANN   size parameters (number of units in the hidden layer) for ANN used for tuning (default: c(2,4,6,8)).  Will be optimised using the method specified in ctrl.ANN (if not available trControl).
##' @param decay.tune.ANN  weight decay parameters used for tuning for ANN (default: c(0.001, 0.01, 0.05, 0.1))  Will be optimised by method specified in ctrl.ANN (if not available trControl).
##' @param maxit.ANN       maximum number of iterations for ANN (default 500) 
##' @param MaxNWts.ANN     The maximum allowable number of weights for ANN (default 10 * (ncol(myBiomodData'at'data.env.var) + 1) + 10 + 1). 
##' @param parallel.ME     logical. If TRUE, the parallel computing is enabled for MAXENT.Phillips 
##' @param numCores.ME     number of cores used to train MAXENT.Phillips 
##' @param Yweights        response points weights. This argument will only affect models that allow case weights. 
##' 
##' @return
##' BIOMOD.models.options object with optimized parameters
##' 
##' @author Frank Breiner \email{frank.breiner@wsl.ch}
##' 
##' @references 
##' Kuhn, Max. 2008. Building predictive models in R using the caret package. \emph{Journal of Statistical Software} \bold{28}, 1-26.
##' Kuhn, Max, and Kjell Johnson. 2013. Applied predictive modeling. New York: Springer.
##' Muscarella, Robert, et al. 2014. ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, \bold{5}, 1198-1205.
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_ModelingOptions}}, \code{\link[caret]{train}}, \code{\link[ENMeval]{ENMevaluate}}, 
##' 
##' @examples
##' \dontrun{
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"))
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio7.grd", 
##'                              package="biomod2"),  
##'                 system.file( "external/bioclim/current/bio11.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio12.grd", 
##'                              package="biomod2"))
##' # 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' ### Duration for turing all models sequential with default settings 
##' ### on 3.4 GHz processor: approx. 45 min tuning all models in parallel
##' ### (on 8 cores) using foreach loops runs much faster: approx. 14 min
##' 
##' #library(doParallel);cl<-makeCluster(8);doParallel::registerDoParallel(cl) 
##' 
##' 
##' time.seq<-system.time(Biomod.tuning <- BIOMOD_tuning(myBiomodData,
##'                                                              env.ME = myExpl,
##'                                                              n.bg.ME = ncell(myExpl)))
##' #stopCluster(cl)
##' 
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##'                                      models = c('RF','CTA'), 
##'                                      models.options = Biomod.tuning$models.options, 
##'                                      NbRunEval=1, 
##'                                      DataSplit=100, 
##'                                      VarImport=0, 
##'                                      models.eval.meth = c('ROC'),
##'                                      do.full.models=FALSE,
##'                                      modeling.id="test")
##' 
##' 
##' #  eval.plot(Biomod.tuning$tune.MAXENT.Phillips at results)
##' par(mfrow=c(1,3))
##' plot(Biomod.tuning$tune.CTA.rpart)
##' plot(Biomod.tuning$tune.CTA.rpart2)
##' plot(Biomod.tuning$tune.RF)
##' }
BIOMOD_tuning <- function(data,
                          models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips'),
                          models.options = BIOMOD_ModelingOptions(),
                          method.ANN = 'avNNet',
                          method.RF = 'rf',
                          method.MARS = 'earth',
                          method.GAM = 'gam',
                          method.GLM = 'glmStepAIC',
                          trControl = NULL,
                          metric = 'ROC',
                          ctrl.CTA = NULL,
                          ctrl.RF = NULL,
                          ctrl.ANN = NULL,
                          ctrl.MARS = NULL,                                  
                          ctrl.FDA = NULL,
                          ctrl.GAM = NULL,
                          ctrl.GBM = NULL,
                          ctrl.GLM = NULL,
                          tuneLength = 30,
                          decay.tune.ANN = c(0.001, 0.01, 0.05, 0.1),
                          size.tune.ANN = c(2,4,6,8),
                          maxit.ANN = 500,
                          MaxNWts.ANN = 10 * (ncol(data@data.env.var) + 1) + 10 + 1,
                          type.GLM = c('simple','quadratic','polynomial','s_smoother'),
                          interaction.GLM = c(0,1),
                          cvmethod.ME = 'randomkfold',
                          overlap.ME = FALSE,
                          kfolds.ME = 10,
                          n.bg.ME = 10000,
                          env.ME = NULL,
                          metric.ME = 'ROC',
                          clamp.ME = TRUE,
                          parallel.ME = FALSE,
                          numCores.ME = NULL,
                          Yweights = NULL){
  
  ## MAXENT: http://cran.r-project.org/web/packages/ENMeval/ENMeval.pdf --> ENMevaluate()
  ## or:    http://cran.r-project.org/web/packages/maxent/maxent.pdf -->  tune.maxent()
  #packages <- NULL
  if(sum(c('GLM','GBM','GAM','CTA','ANN','FDA','MARS','RF','MAXENT.Phillips','SRE') %in% models)>0){if(!isNamespaceLoaded("caret")){requireNamespace("caret", quietly = TRUE)}; 
    if(!isNamespaceLoaded('dplyr')){requireNamespace("dplyr", quietly = TRUE)}; 
    if(is.null(trControl)){trControl <- caret::trainControl(method="cv",summaryFunction=caret::twoClassSummary,classProbs=T, returnData = F)}}
  if("MAXENT.Phillips" %in% models){if(!isNamespaceLoaded('ENMeval')){requireNamespace("ENMeval", quietly = TRUE)}}#;packages<-c(packages,"ENMeval")}  
  # if("MAXENT.Tsuruoka" %in% models){if(!isNamespaceLoaded('maxent')){requireNamespace("maxent", quietly = TRUE)}}#;packages<-c(packages,"maxent")}  
  
  tune.SRE <- tune.GLM <- tune.MAXENT.Phillips <- tune.GAM <- tune.GBM <- tune.CTA.rpart <- tune.CTA.rpart2 <- tune.RF <- tune.ANN <- tune.MARS <- tune.FDA <- NULL
    # tune.MAXENT.Tsuruoka <- NULL
    
  
  resp <- data@data.species
  
  if('SRE' %in% models){
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n",
              "Start tuning SRE\n"))
    #if(trControl$method=="cv"){
    tune.SRE <- NULL
    for (rep in 1:trControl$repeats){
      fold <- dismo::kfold(resp, by = resp, 
                           k = trControl$number)
      for(quant in c(0,0.0125,0.025,0.05,0.1)){
        for (i in 1:trControl$number) {
          DATA <- cbind(1:sum(fold==i),
                        resp[fold==i],
                        sre(Response = resp[fold!=i], 
                            Explanatory = data@data.env.var[fold!=i,], 
                            NewData = data@data.env.var[fold==i,], 
                            Quant=quant, 
                            return_extremcond = FALSE))
          tune.SRE <- rbind(tune.SRE,cbind(
            presence.absence.accuracy(DATA, 
                                      threshold=as.vector(PresenceAbsence::optimal.thresholds(DATA,opt.methods=3)[2],mode="numeric")),quant))
        }
      }}
    t<-aggregate(tune.SRE,by=list(quant = tune.SRE$quant),mean)
    if(metric == 'ROC'){models.options@SRE$quant<-t[which.max(t$AUC),"quant"]} 
    if(metric == 'TSS'){models.options@SRE$quant<-t[which.max(t$sensitivity+t$specificity-1),"quant"]} 
    cat(paste("Finished tuning SRE\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }#}
  
  if(metric == 'ROC' | metric == 'TSS'){resp <- as.factor(ifelse(resp == 1 & !is.na(resp), "Presence", "Absence"))}
  
  if('GBM' %in% models){  
    
    if(is.null(ctrl.GBM)){ctrl.GBM <- trControl}
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n",
              "Start tuning GBM. Start coarse tuning\n"))
    
    tune.grid <- expand.grid(.interaction.depth = seq(2, 8, by = 3),
                             .n.trees = c(500, 1000, 2500),
                             .shrinkage = c(0.001, 0.01, 0.1),
                             .n.minobsinnode = 10)
    
    try(tune.GBM <- caret::train(data@data.env.var, resp,
                                 method = "gbm",
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.GBM,
                                 verbose = FALSE,
                                 weights = Yweights))
    cat("Best optimization of coarse tuning:\n")
    cat(paste(tune.GBM$bestTune,"\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.GBM)){
      cat("Start fine tuning\n")
      
      if(tune.GBM$bestTune$n.trees==2500){
        cat("Best optimization with large trees! Tuning GBM will take a while.\n")
        n.trees <- seq(2500, 10000, by = 2500)}
      
      if(tune.GBM$bestTune$n.trees==1000){n.trees <- seq(750, 2000, by = 250)}
      
      if(tune.GBM$bestTune$n.trees==500){n.trees <- seq(100, 1000, by = 50)}
      
      tune.grid <- expand.grid(.interaction.depth = c(tune.GBM$bestTune$interaction.depth-1,tune.GBM$bestTune$interaction.depth,tune.GBM$bestTune$interaction.depth+1),
                               .n.trees = n.trees,
                               .shrinkage = c(tune.GBM$bestTune$shrinkage/2,tune.GBM$bestTune$shrinkage,tune.GBM$bestTune$shrinkage*5),
                               .n.minobsinnode = 10)
      tune.GBM <- NULL
      try(tune.GBM <- caret::train(data@data.env.var, resp,
                                   method = "gbm",
                                   tuneGrid = tune.grid,
                                   trControl = ctrl.GBM,
                                   verbose = FALSE,
                                   weights = Yweights))  
    }
    cat(paste("\n Finished tuning GBM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.GBM)){
      if(metric == 'TSS'){
        models.options@GBM$n.trees <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"n.trees"]
        models.options@GBM$interaction.depth <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"interaction.depth"]    
        models.options@GBM$shrinkage <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"shrinkage"] 
      }else{
        models.options@GBM$n.trees <- tune.GBM$bestTune$n.trees
        models.options@GBM$interaction.depth <- tune.GBM$bestTune$interaction.depth    
        models.options@GBM$shrinkage <- tune.GBM$bestTune$shrinkage 
      }}else{ if('GBM' %in% models){cat("Tuning GBM failed!"); tune.GBM <- "FAILED"}}
    
  }
  
  if('RF' %in% models){
    cat("Start tuning RF\n")
    
    if(is.null(ctrl.RF)){ctrl.RF <- trControl}
    tuneLength.rf <- min(tuneLength,ncol(data@data.env.var))
    
    ## give both mtry as bestTune
    try(tune.RF <- caret::train(data@data.env.var, resp,
                                method = method.RF,
                                tuneLength = tuneLength.rf,
                                trControl = ctrl.RF,
                                metric = metric,
                                weights = Yweights))
    cat(paste("Finished tuning RF\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    if(!is.null(tune.RF)){
      if(metric == 'TSS'){
        models.options@RF$mtry <- tune.RF$results[which.max(apply(tune.RF$results[,c("Sens","Spec")],1,sum)-1),"mtry"]
      }else{
        models.options@RF$mtry <- tune.RF$bestTune$mtry
      }}else{ if('RF' %in% models){cat("Tuning RF failed!"); tune.RF <- "FAILED"}}
  }
  
  if('ANN' %in% models){
    cat("Start tuning ANN\n")
    
    if(is.null(ctrl.ANN)){ctrl.ANN <- trControl}
    ## already tuning: 
    # size: optimised by cross validation based on model AUC (NbCv cross validation; tested size will be the following c(2,4,6, 8))
    # decay: optimised by cross validation on model AUC (NbCv cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1)).
    # could increase maxit from 200 to 500
    # a nice option would be to use model averaging for ann: avNNet in package(caret)
    
    ## Create a specific candidate set of models to evaluate:
    tune.grid <- expand.grid(.decay = decay.tune.ANN,
                             .size = size.tune.ANN,
                             .bag = FALSE)
    
    try(tune.ANN <- caret::train(data@data.env.var, resp, 
                                 method = method.ANN,
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.ANN,
                                 ## Automatically standardize data prior to modeling
                                 ## and prediction
                                 preProc = c("center", "scale"),
                                 linout = TRUE,
                                 trace = FALSE,
                                 MaxNWts.ANN = MaxNWts.ANN,
                                 maxit = maxit.ANN,
                                 metric = metric,
                                 weights = Yweights))
    cat(paste("Finished tuning ANN\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    if(!is.null(tune.ANN)){
      if(metric == 'TSS'){
        models.options@ANN$size <- tune.ANN$results[which.max(apply(tune.ANN$results[,c("Sens","Spec")],1,sum)-1),"size"]    
        models.options@ANN$decay <- tune.ANN$results[which.max(apply(tune.ANN$results[,c("Sens","Spec")],1,sum)-1),"decay"]    
        models.options@ANN$maxit <- maxit.ANN
      }else{
        models.options@ANN$size <- tune.ANN$bestTune$size
        models.options@ANN$decay <- tune.ANN$bestTune$decay 
        models.options@ANN$maxit <- maxit.ANN
      }}else{ if('ANN' %in% models){cat("Tuning ANN failed!"); tune.ANN <- "FAILED"}}
    
  }
  
  if('GAM' %in% models){
    cat("Start tuning GAM\n")
    
    if(is.null(ctrl.GAM)){ctrl.GAM <- trControl}
    
    try(tune.GAM <-   caret::train(data@data.env.var, resp, 
                                   method = method.GAM,
                                   trControl = ctrl.GAM,
                                   weights = Yweights))
    cat(paste("Finished tuning GAM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.GAM)){
      if(metric == 'TSS'){
        models.options@GAM$select <- tune.GAM$results[which.max(apply(tune.GAM$results[,c("Sens","Spec")],1,sum)-1),"select"]    
        models.options@GAM$method <- as.character(tune.GAM$results[which.max(apply(tune.GAM$results[,c("Sens","Spec")],1,sum)-1),"method"])    
      }else {
        models.options@GAM$select <- tune.GAM$bestTune$select
        models.options@GAM$method <- as.character(tune.GAM$bestTune$method)
      }}else{ if('GAM' %in% models){cat("Tuning GAM failed!"); tune.GAM <- "FAILED"}}
  }
  
  if('MARS' %in% models){
    cat("Start tuning MARS\n")
    
    if(is.null(ctrl.MARS)){ctrl.MARS <- trControl}
    
    if(is.null(models.options@MARS$nk)){nprune <- 2:max(21, 2 * ncol(data@data.env.var) + 1)
    }else{
      nprune <- 2:min(models.options@MARS$nk,38)}
    tune.grid <- expand.grid(.degree = 1:2, .nprune = nprune)
    try(tune.MARS <-   caret::train(data@data.env.var, resp, 
                                    method = method.MARS,
                                    tuneGrid = tune.grid,
                                    trControl = ctrl.MARS,
                                    weights = Yweights))
    
    cat(paste("Finished tuning MARS\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.MARS)){
      if(metric == 'TSS'){
        if("degree" %in% names(models.options@MARS)){
          models.options@MARS$degree <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"degree"]    
        }else{
          models.options@MARS$interaction.level <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"degree"]-1    
        }
        models.options@MARS$nprune <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"nprune"]    
      }else{
        if("degree" %in% names(models.options@MARS)){
          models.options@MARS$degree <- tune.MARS$bestTune$degree
        }else{
          models.options@MARS$interaction.level <- tune.MARS$bestTune$degree-1    
        }
        models.options@MARS$nprune <- tune.MARS$bestTune$nprune
      }}else{ if('MARS' %in% models){cat("Tuning MARS failed!"); tune.MARS <- "FAILED"}}
  }
  
  
  if('GLM' %in% models){
    cat("Start tuning GLM\n")
    
    if(is.null(ctrl.GLM)){ctrl.GLM <- trControl}
    if("s_smoother" %in% type.GLM){requireNamespace("gam", quietly = TRUE)}
    
    fm<-list()
    GLM.results<-NULL 
    i<-0
    for(type in type.GLM){
      for(IA in interaction.GLM){
        i<-i+1
        try(tune.GLM <-   caret::train( makeFormula("resp",data@data.env.var, type= type,interaction.level = IA),
                                        data=cbind(data@data.env.var,resp=resp),
                                        method = method.GLM,
                                        trControl = ctrl.GLM,
                                        weights = Yweights))  
        try(GLM.results <-  rbind(GLM.results,cbind(tune.GLM$results,il=IA,type=type)))
        try(fm[[i]] <- formula(tune.GLM$finalModel))
      }
    } 
    
    glm.best<-which.max(GLM.results$ROC)
    models.options@GLM$interaction.level <- GLM.results[glm.best,"il"]     
    models.options@GLM$type <- as.character(GLM.results[glm.best,"type"])
    models.options@GLM$myFormula <- formula(paste(data@sp.name,"~",gsub("`","",as.character(fm[[glm.best]])[3])))
    models.options@GLM$test <- "none" 
    
    cat(paste("Finished tuning GLM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }      
  
  
  if('FDA' %in% models){
 
    cat("Start tuning FDA\n")
    
    if(is.null(ctrl.FDA)){ctrl.FDA <- trControl}
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = 2:38)
    try(tune.FDA <- caret::train(data@data.env.var, factor(resp), 
                                 method = "fda",
                                 tuneGrid = tune.grid,                  
                                 trControl = ctrl.FDA,
                                 weights = Yweights))
    cat(paste("Finished tuning FDA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.FDA)){
      #models.options@FDA$method <- "earth"   
      if(metric == 'TSS'){
        models.options@FDA$add_args <- list(degree=tune.FDA$results[which.max(apply(tune.FDA$results[,c("Sens","Spec")],1,sum)-1),"degree"],   
                                            nprune=tune.FDA$results[which.max(apply(tune.FDA$results[,c("Sens","Spec")],1,sum)-1),"nprune"])    
      }else{
        models.options@FDA$add_args <- list(degree=tune.FDA$bestTune$degree,nprune=tune.FDA$bestTune$nprune)
      }}else{ if('FDA' %in% models){cat("Tuning FDA failed!"); tune.FDA <- "FAILED"}}
  }
  
  if('CTA' %in% models){
    cat("Start tuning CTA\n")
    
    if(is.null(ctrl.CTA)){ctrl.CTA <- trControl}    
    
    cat("Tuning Complexity Parameter")    
    try(tune.CTA.rpart <- caret::train(data@data.env.var, resp, 
                                       method = "rpart",
                                       tuneLength = tuneLength,
                                       trControl = ctrl.CTA,
                                       metric=metric,
                                       weights = Yweights))
    
    cat("Tuning Max Tree Depth")
    try(tune.CTA.rpart2 <-  caret::train(data@data.env.var, resp,
                                         method = "rpart2",
                                         tuneLength = tuneLength,
                                         trControl = ctrl.CTA,
                                         metric=metric,
                                         weights = Yweights))
    cat(paste("Finished tuning CTA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.CTA.rpart)){
      if(metric == 'TSS'){
        models.options@CTA$control$cp <- tune.CTA.rpart$results[which.max(apply(tune.CTA.rpart$results[,c("Sens","Spec")],1,sum)-1),"cp"]    
      }else{
        models.options@CTA$control$cp <- tune.CTA.rpart$bestTune      
      }}else{ if('CTA' %in% models){cat("Tuning CTA cp failed!"); tune.CTA.rpart <- "FAILED"}}
    
    if(!is.null(tune.CTA.rpart2)){
      if(metric == 'TSS'){
        models.options@CTA$control$maxdepth <- tune.CTA.rpart2$results[which.max(apply(tune.CTA.rpart2$results[,c("Sens","Spec")],1,sum)-1),"maxdepth"]    
      }else{
        models.options@CTA$control$maxdepth <- tune.CTA.rpart2$bestTune      
      }}else{ if('CTA' %in% models){cat("Tuning CTA maxdepth failed!"); tune.CTA.rpart2 <- "FAILED"}}
  }
  
  if('MAXENT.Phillips' %in% models){
    cat("Start tuning MAXENT.Phillips\n")
    if(cvmethod.ME != 'randomkfold'){kfolds.ME <- NA}
    try(tune.MAXENT.Phillips <- tuning.maxent(pres=data@data.env.var[data@data.species==1 & !is.na(data@data.species),],
                                              bg= data@data.env.var[data@data.species==0 | is.na(data@data.species),],
                                              method=cvmethod.ME, kfolds = kfolds.ME,#env.ME,
                                              bin.output=TRUE, clamp=clamp.ME, parallel = parallel.ME, numCores = numCores.ME,
                                              categoricals=NULL))
    cat(paste("Finished tuning MAXENT.Phillips\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.MAXENT.Phillips)){
      if(metric.ME=="ROC"){metric.ME <- "Mean.AUC"}
      if(!metric.ME %in% c("Mean.AUC", "Mean.AUC.DIFF", "Mean.ORmin", "Mean.OR10", "AICc")){metric.ME <- "Mean.AUC"; cat("Invalid metric.ME argument! metric.ME was set to Mean.AUC")}
      if(metric.ME == 'Mean.AUC'){
        models.options@MAXENT.Phillips$linear <- grepl("L",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$quadratic <- grepl("Q",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$hinge <- grepl("H",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$product <- grepl("P",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$threshold <- grepl("T",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$betamultiplier <- tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"rm"]  
      }else {       
        models.options@MAXENT.Phillips$linear <- grepl("L",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$quadratic <- grepl("Q",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$hinge <- grepl("H",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$product <- grepl("P",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$threshold <- grepl("T",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$betamultiplier <- tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"rm"]   
      }}else{ if('MAXENT.Phillips' %in% models){cat("Tuning MAXENT.Phillips failed!"); tune.MAXENT.Phillips <- "FAILED"}}
  }
  
  
  # if('MAXENT.Tsuruoka' %in% models){
  #   cat("Start tuning MAXENT.Tsuruoka\n")
  #   try(tune.MAXENT.Tsuruoka <- as.data.frame(tune.maxent(data@data.env.var,data@data.species,nfold=kfolds.ME,showall=T)))
  #   cat(paste("Finished tuning MAXENT.Tsuruoka\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  #   
  #   if(!is.null(tune.MAXENT.Tsuruoka)){
  #     models.options@MAXENT.Tsuruoka$l1_regularizer <- tune.MAXENT.Tsuruoka$l1_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     models.options@MAXENT.Tsuruoka$l2_regularizer <- tune.MAXENT.Tsuruoka$l2_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     models.options@MAXENT.Tsuruoka$use_sgd <- ifelse(tune.MAXENT.Tsuruoka[which.max(tune.MAXENT.Tsuruoka$accuracy),]$use_sgd==0,F,T)
  #     models.options@MAXENT.Tsuruoka$set_heldout <- tune.MAXENT.Tsuruoka$set_heldout[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #   } else { if('MAXENT.Tsuruoka' %in% models){cat("Tuning MAXENT.Tsuruoka failed!"); tune.MAXENT.Tsuruoka <- "FAILED"}}
  # }  
  
  
  return(list(models.options=models.options, tune.SRE =tune.SRE,  tune.CTA.rpart = tune.CTA.rpart, tune.CTA.rpart2 = tune.CTA.rpart2,
              tune.RF = tune.RF, tune.ANN = tune.ANN,  tune.MARS = tune.MARS, tune.FDA = tune.FDA, tune.GBM=tune.GBM,
              tune.GAM = tune.GAM, tune.MAXENT.Phillips = tune.MAXENT.Phillips, 
              # tune.MAXENT.Tsuruoka = tune.MAXENT.Tsuruoka, 
              tune.GLM=tune.GLM))
}


##
#### Modified tuning function from the ENMeval package to tune MAXENT.Phillips (internal function for BIOMOD_tuning)

tuning.maxent <-
  function (occ, env=NULL,pres=NULL, bg=NULL, bg.coords=NULL, occ.grp=NULL, bg.grp=NULL, method=NULL, maxent.args, 
            args.lab, categoricals=NULL, aggregation.factor=c(2,2), kfolds=NA, bin.output=FALSE, 
            clamp, rasterPreds=FALSE, parallel=FALSE, numCores=NULL,RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")) 
  {
    ###NEW
    requireNamespace("ENMeval", quietly = TRUE)
    maxent.args <- ENMeval::make.args(RMvalues, fc)
    args.lab <- ENMeval::make.args(RMvalues, fc, labels = TRUE)
    if(!is.null(pres)){occ<-pres}
    if(!is.null(bg)){bg.coords<-bg}
    #####
    noccs <- nrow(occ)
    if (method == "checkerboard1") 
      group.data <- ENMeval::get.checkerboard1(occ, env, bg.coords, 
                                               aggregation.factor)
    if (method == "checkerboard2") 
      group.data <- ENMeval::get.checkerboard2(occ, env, bg.coords, 
                                               aggregation.factor)
    if (method == "block") 
      group.data <- ENMeval::get.block(occ, bg.coords)
    if (method == "jackknife") 
      group.data <- ENMeval::get.jackknife(occ, bg.coords)
    if (method == "randomkfold") 
      group.data <- ENMeval::get.randomkfold(occ, bg.coords, kfolds)
    if (method == "user") 
      group.data <- ENMeval::get.user(occ.grp, bg.grp)
    nk <- length(unique(group.data$occ.grp))
    ###NEW
    if(is.null(pres)){
      pres <- as.data.frame(extract(env, occ))}
    if(is.null(bg)){
      bg <- as.data.frame(extract(env, bg.coords))}
    #####
    if (any(is.na(pres))) {
      message("Warning: some predictors variables are NA at some occurrence points")
    }
    if (any(is.na((bg)))) {
      message("Warning: some predictors variables are NA at some background points")
    }
    if (!is.null(categoricals)) {
      for (i in 1:length(categoricals)) {
        pres[, categoricals[i]] <- as.factor(pres[, categoricals[i]])
        bg[, categoricals[i]] <- as.factor(bg[, categoricals[i]])
      }
    }
    tune <- function() {
      if (length(maxent.args) > 1 & !parallel) {
        setTxtProgressBar(pb, i)
      }
      x <- rbind(pres, bg)
      p <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
      tmpfolder <- tempfile()
      full.mod <- dismo::maxent(x, p, args = maxent.args[[i]], factors = categoricals, 
                                path = tmpfolder)
      pred.args <- c("outputformat=raw", ifelse(clamp == TRUE, 
                                                "doclamp=true", "doclamp=false"))
      if (rasterPreds == TRUE) {
        predictive.map <- predict(full.mod, env, args = pred.args)
      }
      else {
        predictive.map <- stack()
      }
      AUC.TEST <- double()
      AUC.DIFF <- double()
      OR10 <- double()
      ORmin <- double()
      for (k in 1:nk) {
        train.val <- pres[group.data$occ.grp != k, ]
        test.val <- pres[group.data$occ.grp == k, ]
        bg.val <- bg[group.data$bg.grp != k, ]
        x <- rbind(train.val, bg.val)
        p <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))
        mod <- dismo::maxent(x, p, args = maxent.args[[i]], factors = categoricals, 
                             path = tmpfolder)
        ### Specify dismo!!! Problems with biomod!
        AUC.TEST[k] <- dismo::evaluate(test.val, bg, mod)@auc
        AUC.DIFF[k] <- max(0, dismo::evaluate(train.val, bg, mod)@auc - 
                             AUC.TEST[k])
        ########
        p.train <- predict(mod, train.val, args = pred.args)
        p.test <- predict(mod, test.val, args = pred.args)
        if (nrow(train.val) < 10) {
          n90 <- floor(nrow(train.val) * 0.9)
        }
        else {
          n90 <- ceiling(nrow(train.val) * 0.9)
        }
        train.thr.10 <- rev(sort(p.train))[n90]
        OR10[k] <- mean(p.test < train.thr.10)
        train.thr.min <- min(p.train)
        ORmin[k] <- mean(p.test < train.thr.min)
      }
      unlink(tmpfolder, recursive = TRUE)
      stats <- c(AUC.DIFF, AUC.TEST, OR10, ORmin)
      return(list(full.mod, stats, predictive.map))
    }
    if (parallel == TRUE) {
      requireNamespace("foreach")
      allCores <- detectCores()
      if (is.null(numCores)) {
        numCores <- allCores
      }
      c1 <- makeCluster(numCores)
      doParallel::registerDoParallel(c1)
      numCoresUsed <- foreach::getDoParWorkers()
      message(paste("Of", allCores, "total cores using", numCoresUsed))
      message("Running in parallel...")
      out <- foreach::foreach(i = seq_len(length(maxent.args)), .packages = c("dismo", 
                                                                              "raster", "ENMeval")) %dopar% {
                                                                                tune()
                                                                              }
      stopCluster(c1)
    }
    else {
      pb <- txtProgressBar(0, length(maxent.args), style = 3)
      out <- list()
      for (i in 1:length(maxent.args)) {
        out[[i]] <- tune()
      }
      close(pb)
    }
    full.mods <- sapply(out, function(x) x[[1]])
    statsTbl <- as.data.frame(t(sapply(out, function(x) x[[2]])))
    if (rasterPreds) {
      predictive.maps <- stack(sapply(out, function(x) x[[3]]))
    }
    else {
      predictive.maps <- stack()
    }
    AUC.DIFF <- statsTbl[, 1:nk]
    AUC.TEST <- statsTbl[, (nk + 1):(2 * nk)]
    OR10 <- statsTbl[, ((2 * nk) + 1):(3 * nk)]
    ORmin <- statsTbl[, ((3 * nk) + 1):(4 * nk)]
    names(AUC.DIFF) <- paste("AUC.DIFF_bin", 1:nk, sep = ".")
    Mean.AUC.DIFF <- rowMeans(AUC.DIFF)
    Var.AUC.DIFF <- ENMeval::corrected.var(AUC.DIFF, noccs)
    names(AUC.TEST) <- paste("AUC_bin", 1:nk, sep = ".")
    Mean.AUC <- rowMeans(AUC.TEST)
    Var.AUC <- ENMeval::corrected.var(AUC.TEST, noccs)
    names(OR10) <- paste("OR10_bin", 1:nk, sep = ".")
    Mean.OR10 <- rowMeans(OR10)
    Var.OR10 <- apply(OR10, 1, var)
    names(ORmin) <- paste("ORmin_bin", 1:nk, sep = ".")
    Mean.ORmin <- rowMeans(ORmin)
    Var.ORmin <- apply(ORmin, 1, var)
    full.AUC <- double()
    for (i in 1:length(full.mods)) full.AUC[i] <- full.mods[[i]]@results[5]
    nparm <- numeric()
    for (i in 1:length(full.mods)) nparm[i] <- ENMeval::get.params(full.mods[[i]])
    if (rasterPreds == TRUE) {
      aicc <- ENMeval::calc.aicc(nparm, occ, predictive.maps)
    }
    else {
      aicc <- rep(NaN, length(full.AUC))
    }
    features <- args.lab[[1]]
    rm <- args.lab[[2]]
    settings <- paste(args.lab[[1]], args.lab[[2]], sep = "_")
    res <- data.frame(settings, features, rm, full.AUC, Mean.AUC, 
                      Var.AUC, Mean.AUC.DIFF, Var.AUC.DIFF, Mean.OR10, Var.OR10, 
                      Mean.ORmin, Var.ORmin, nparm, aicc)
    if (bin.output == TRUE) {
      res <- as.data.frame(cbind(res, AUC.TEST, AUC.DIFF, OR10, 
                                 ORmin))
    }
    if (rasterPreds == TRUE) {
      names(predictive.maps) <- settings
    }
    results <- ENMeval::ENMevaluation(results = res, predictions = predictive.maps, 
                                      models = full.mods, partition.method = method, occ.pts = occ, 
                                      occ.grp = group.data[[1]], bg.pts = bg.coords, bg.grp = group.data[[2]])
    return(results)
  }



.CleverCut <- function(x){
#### old version
#   switch(EXPR=x,
#          '1' = return(c(1,1)),
#          '2' = return(c(1,2)),
#          '3' = return(c(2,2)),
#          '4' = return(c(2,2)),
#          '5' = return(c(2,3)),
#          '6' = return(c(2,3)),
#          '7' = return(c(3,3)),
#          '8' = return(c(3,3)),
#          return(c(3,3)))
  
  nb_col = ceiling(sqrt(x))
  nb_row = ceiling(x/nb_col)
  return(c(nb_row,nb_col))
}

.bmCat <- function(x=NULL,...){
  if(is.null(x)){
    cat("\n")
    cat(paste(rep("-=", round(.Options$width/2) ), collapse=""))
    cat("\n")
  } else{
    x.length = nchar(x) + 2
    y.length = (.Options$width - x.length) / 2
    cat("\n")
    cat(paste(rep("-=", round(y.length/2) ), collapse=""), x, paste(rep("-=", round(y.length/2) ), collapse=""))
    cat("\n")
  }
}
'CustomIndexMaker' <- function(){
  # 1. test if modif has ever been done
#   return(FALSE)

  if(file.exists(paste(system.file("",package='biomod2'),
                 .Platform$file.sep, "HasBeenCustom.txt", sep=""))){
    return(FALSE)
  } else{
    cat("\nbiomod2 first load...")

    if(file.exists(file.path(system.file("doc",package='biomod2'), "html","00Index.html"))){
      cat("\ncustom help index setting up...")
      # get currentindex files
      old.index <- system.file("html/00Index.html",package='biomod2')

      # get customed index file
      new.index <- system.file("doc/html/00Index.html",package='biomod2')

      # update version number
      pkg.version <- read.dcf(file=system.file("DESCRIPTION", package='biomod2'),
               fields="Version")

      indexTmp <- readLines(new.index)
      indexTmp[grepl("version ",indexTmp)] <- paste("version ", pkg.version, sep="")
      cat(indexTmp, sep="\n", file=new.index, append=FALSE)

      # and replace it
      file.copy(from=old.index,to=paste(tools::file_path_sans_ext(old.index),"_DEFAULT.html", sep=""))
      file.copy(from=new.index,to=old.index,overwrite=TRUE)
    }

#     if(!file.exists(file.path(system.file("",package='biomod2'), "Meta", "vignette.rds"))){
#       cat("\ncustom vignettes files setting up...")
#       # get vignettes files
#       vignettes_R <- list.files(path=system.file("doc",package='biomod2'),pattern=".R",ignore.case=TRUE)
#       vignettes_PDF <- list.files(path=system.file("doc",package='biomod2'),pattern=".pdf",ignore.case=TRUE)
#
#       # create the .Rnw & vignette.rds file
#       vignettes_Rnw <- union(gsub(".R","",vignettes_R), gsub(".pdf","",vignettes_PDF))
#       vignettes_Rnw <- paste(vignettes_Rnw,".Rnw", sep="")
#       invisible(file.create(file.path(system.file("doc",package='biomod2'), vignettes_Rnw), showWarnings=F))
#
#       vignettes_rds <- file.path(system.file("doc",package='biomod2'),"html","vignettes_rds.csv")
#       if(file.exists(vignettes_rds)){
#         vignettes_rds <- read.csv(vignettes_rds, colClasses = "character")
#         saveRDS(vignettes_rds,file=file.path(system.file("Meta",package='biomod2'),"vignette.rds"))
#       }
#     }
    # create a file that will indicate that package customing has ever been done
    invisible(file.create(file.path(system.file("",package='biomod2'), "HasBeenCustom.txt"), showWarnings=F))
    return(TRUE)
  }
}




.onAttach <- function(libname, pkgname) {
  if(file.exists(system.file("DESCRIPTION", package=pkgname))){
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")

    if(CustomIndexMaker()){
      packageStartupMessage("Custom index built!")
    }

    packageStartupMessage(paste(pkgname, RFver, "loaded.\n\nType browseVignettes(package='biomod2') to access directly biomod2 vignettes."))
  }
}
.CV.nnet = function(Input, Target, size=c(2,4,6, 8), decay=c(0.001, 0.01, 0.05, 0.1), maxit=200, nbCV=5, W=NULL){
#   require(pROC, quietly=T)

  Eval = data.frame(matrix(0, ncol=3, nrow=16, dimnames=list(NULL, c("Size", "Decay", "AUC"))))
  Eval[,1] = rep(size,4)
  Eval[,2] = rep(decay, each=4)
  for(i in 1:nbCV){
      set.seed(555)
      Samp = SampleMat2(Target, 0.5)
      
      if(is.null(W)){
          Eval[,3] = Eval[,3] + apply(Eval[,1:2], 1, Samp, Target, Input, FUN=function(x, Samp, Target, Input){
            nn = nnet(eval(parse(text = paste("Target[Samp$calibration]",
                  paste(.scopeExpSyst(Input[1:10, ,drop=FALSE], "GBM"), collapse = "")))),data=Input[Samp$calibration, ,drop=FALSE],
                  size = x[1], decay = x[2], maxit = maxit, trace = FALSE)
            AUC = as.numeric(pROC::auc(pROC::roc(Target[Samp$evaluation], predict(nn, Input[Samp$evaluation,,drop=FALSE]))))
            return(AUC)
          })
      } else{
          Eval[,3] = 
            Eval[,3] + 
            apply(
              Eval[,1:2], 1, Samp, Target, Input, W, 
              FUN = function(x, Samp, Target, Input, W){
                nn = 
                  nnet(
                    eval(parse(text = paste("Target[Samp$calibration]", paste(.scopeExpSyst(Input[1:10, ,drop=FALSE], "GBM"), collapse = "")))),
                    data = Input[Samp$calibration, ,drop=FALSE],
                    weights = W[Samp$calibration], 
                    size = x[1], 
                    decay = x[2], 
                    maxit = maxit, 
                    trace = FALSE
                  )
                AUC <- 
                  as.numeric(
                    pROC::auc(
                      pROC::roc(
                        Target[Samp$evaluation], 
                        as.numeric(predict(nn, Input[Samp$evaluation,,drop=FALSE])),
                        levels = c(0, 1),
                        direction = '<'
                      )
                    )
                  )
            return(AUC)
          })
      }
  }
  Eval[, 3] = Eval[, 3] / nbCV
  z = which.max(Eval[, 3])
  return(Eval[z, 1:2])
}
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Compilation of deprecated function that will be removed from
# the package one day or another
# Damien G. - 26/06/2013
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setGeneric("getModelsPrediction",
           function(obj,...){
             standardGeneric("getModelsPrediction")
           })

setMethod("getModelsPrediction", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_predictions(obj, eval_data=FALSE)")
            if(!as.data.frame){
              if(obj@models.prediction@inMemory ){
                return(obj@models.prediction@val)
              } else{
                if(obj@models.prediction@link != ''){
                  #                   load(obj@models.prediction@link)
                  #                   return(models.prediction)

                  return(get(load(obj@models.prediction@link)))
                } else{ return(NULL) }
              }
            } else {
              if(obj@models.prediction@inMemory ){
                mod.pred <- as.data.frame(obj@models.prediction@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE),
                                                 function(x){
                                                   return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                 }))
                return(mod.pred)
              } else{
                if(obj@models.prediction@link != ''){
                  #                   load(obj@models.prediction@link)
                  #                   mod.pred <- as.data.frame(models.prediction)
                  mod.pred <- as.data.frame(get(load(obj@models.prediction@link)))
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE),
                                                   function(x){
                                                     return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
                  return(mod.pred)
                } else{ return(NULL) }
              }

            }
          }
)


setGeneric("getModelsPredictionEval",
           function(obj,...){
             standardGeneric("getModelsPredictionEval")
           })

setMethod("getModelsPredictionEval", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_predictions(obj, eval_data=TRUE)")
            if(!as.data.frame){
              if(obj@models.prediction.eval@inMemory ){
                return(obj@models.prediction.eval@val)
              } else{
                if(obj@models.prediction.eval@link != ''){
                  models.prediction.eval <- get(load(obj@models.prediction.eval@link))
                  return(models.prediction.eval)
                } else{ return(NULL) }
              }
            } else {
              if(obj@models.prediction.eval@inMemory ){
                mod.pred <- as.data.frame(obj@models.prediction.eval@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE),
                                                 function(x){
                                                   return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                 }))
                return(mod.pred)
              } else{
                if(obj@models.prediction.eval@link != ''){
                  load(obj@models.prediction.eval@link)
                  mod.pred <- as.data.frame(models.prediction.eval)
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE),
                                                   function(x){
                                                     return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
                  return(mod.pred)
                } else{ return(NULL) }
              }

            }
          }
)


setGeneric("getModelsEvaluations",
           function(obj,...){
             standardGeneric("getModelsEvaluations")
           })

setMethod("getModelsEvaluations", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_evaluations(obj)")
            return(get_evaluations(obj=obj))
          }
)


setGeneric("getModelsVarImport",
           function(obj,...){
             standardGeneric("getModelsVarImport")
           })

setMethod("getModelsVarImport", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_variables_importance(obj)")
            if(obj@variables.importances@inMemory ){
              return(obj@variables.importances@val)
            } else{
              if(obj@variables.importances@link != ''){
                #                 load(obj@variables.importances@link)
                #                 return(variables.importances)
                return(get(load(obj@variables.importances@link)))
              } else{ return(NA) }
            }
          }
)

setGeneric("getModelsOptions",
           function(obj,...){
             standardGeneric("getModelsOptions")
           })

setMethod("getModelsOptions", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_options(obj)")
            if(obj@models.options@inMemory ){
              return(obj@models.options@val)
            } else{
              if(obj@models.options@link != ''){
                #                 load(obj@models.options@link)
                #                 return(models.options)
                return(get(load(obj@models.options@link)))
              } else{ return(NA) }
            }
          }
)

setGeneric("getModelsInputData",
           function(obj, ...){
             standardGeneric("getModelsInputData")
           })

setMethod("getModelsInputData", "BIOMOD.models.out",
          function(obj, subinfo = NULL){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_predictions(obj)")
            if(is.null(subinfo)){
              if(obj@formated.input.data@inMemory ){
                return(obj@formated.input.data@val)
              } else{
                if(obj@formated.input.data@link != ''){
                  data <- get(load(obj@formated.input.data@link))
                  return(data)
                } else{ cat("\n***"); return(NA) }
              }
            } else if(subinfo == 'MinMax'){
              return(apply(getModelsInputData(obj)@data.env.var,2, function(x){
                if(is.numeric(x)){
                  return( list(min = min(x,na.rm=T), max = max(x, na.rm=T) ) )
                } else if(is.factor(x)){
                  return(list(levels = levels(x)))
                }
              }) )
            } else if(subinfo == 'expl.var'){
              return(as.data.frame(getModelsInputData(obj)@data.env.var))
            } else if(subinfo == 'expl.var.names'){
              return(obj@expl.var.names)
            } else if(subinfo == 'resp.var'){
              return(as.numeric(getModelsInputData(obj)@data.species))
            } else if(subinfo == 'eval.resp.var'){
              return(as.numeric(getModelsInputData(obj)@eval.data.species))
            } else if(subinfo == 'eval.expl.var'){
              return(as.data.frame(getModelsInputData(obj)@eval.data.env.var))
            } else{
              stop("Unknown subinfo tag")
            }

          }
)

setGeneric("getModelsBuiltModels",
           function(obj,...){
             standardGeneric("getModelsBuiltModels")
           })

setMethod("getModelsBuiltModels", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_built_models(obj)")
            return(obj@models.computed)
          }
)

setGeneric("getProjection",
           function(obj, ...){
             standardGeneric("getProjection")
           })

setMethod("getProjection", "BIOMOD.projection.out",
          function(obj, model = NULL, as.data.frame = FALSE){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_predictions(obj)")
            if(!as.data.frame & is.null(model)){
              if(obj@proj@inMemory ){
                return(obj@proj@val)
              } else {
                if( grepl(".RData", obj@proj@link) ){
                  return(get(load(obj@proj@link)))
                } else if(grepl(".grd", obj@proj@link) | grepl(".img", obj@proj@link)){
                  return(raster::stack(obj@proj@link, RAT=FALSE))
                } else {
                  filesToLoad <- list.files(path=sub("/individual_projections","", obj@proj@link), full.names=T)
                  toMatch <- c('.grd$','.img$')
                  filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
                  if(length(filesToLoad)){
                    return(raster::stack(filesToLoad[1], RAT=FALSE))
                  } else {
                    filesToLoad <- list.files(path=obj@proj@link, full.names=T)
                    toMatch <- c('.grd$','.img$')
                    filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
                    toMatch <- obj@models.projected
                    filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
                    proj <- raster::stack(filesToLoad, RAT=FALSE)
                    toMatch <- c(obj@proj@link,".img$",'.grd$', .Platform$file.sep)
                    names(proj) <- gsub(pattern=paste(toMatch,collapse="|"), "", filesToLoad)
                    return(proj)
                  }
                }
              }
            } else if(as.data.frame){
              if(obj@proj@inMemory ){
                proj <- as.data.frame(obj@proj@val)
                names(proj) <- unlist(lapply(strsplit(names(proj),".", fixed=TRUE),
                                             function(x){
                                               return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                             }))
                return(proj)
              } else{
                if(obj@proj@link != ''){
                  load(obj@proj@link)
                  project <- as.data.frame(proj)
                  names(project) <- unlist(lapply(strsplit(names(project),".", fixed=TRUE),
                                                  function(x){
                                                    return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                  }))
                  return(project)
                } else{ return(NA) }
              }
            }

          }
)


setGeneric("getEMalgos",
           function(obj,...){
             standardGeneric("getEMalgos")
           })

setMethod("getEMalgos", "BIOMOD.EnsembleModeling.out",
          function(obj, model){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_built_models(obj)")
            return(get_built_models(obj))

          }
)

setGeneric("getEM_needed_models",
           function(obj, ...){
             standardGeneric("getEM_needed_models")
           })

setMethod("getEM_needed_models", "BIOMOD.EnsembleModeling.out",
          function(obj, subset='all', ...){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_needed_models(obj)")
            add.args <- list(...)
            needed_models <- lapply(obj@em.models, function(x){
              return(x@model)
            })
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          }
)


setGeneric("getEMkeptModels",
           function(obj,...){
             standardGeneric("getEMkeptModels")
           })

setMethod("getEMkeptModels", "BIOMOD.EnsembleModeling.out",
          function(obj, model){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_kept_models(obj)")
            if(is.character(model) | is.numeric(model)){
              return(obj@em.res[[model]]$em.models.kept)
            } else{
              return(NULL)
            }

          }
)

setGeneric("getEMeval",
           function(obj, ...){
             standardGeneric("getEMeval")
           })

setMethod("getEMeval", "BIOMOD.EnsembleModeling.out",
          function(obj, model=NULL, met=NULL){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_evaluations(obj)")
            return(get_evaluations(obj=obj, model=model, met=met))


          }
)

setGeneric("getEMbuiltModels",
           function(obj,...){
             standardGeneric("getEMbuiltModels")
           })

setMethod("getEMbuiltModels", "BIOMOD.EnsembleModeling.out",
          function(obj){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_built_models(obj)")
            return(obj@em.computed)
          })

setGeneric( "getFormalModel",
            def = function(obj,...){
              standardGeneric( "getFormalModel" )
            } )

setGeneric( "getScalingModel",
            def = function(obj,...){
              standardGeneric( "getScalingModel" )
            } )

setMethod('getFormalModel', signature('biomod2_model'),
          function(obj){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_formal_model(obj)")
            return(obj@model)
          })

setMethod('getScalingModel', signature('biomod2_model'),
          function(obj){
            cat("\n!! deprecated function that will be removed in next package update")
            cat("\n please use get_scaling_model(obj)")
            return(obj@scaling_model)
          })


evaluate <- function(model, data, stat, as.array=FALSE){
  ## output initialisation
  eval <- NULL

  if(inherits(model, "BIOMOD.models.out") | inherits(model,"BIOMOD.EnsembleModeling.out")){
    eval <- .evaluate.biomod2.models.out(mod=model,data=data,stat=stat)
    if(as.array) eval <- LIST_to_ARRAY(eval)
  } else if(inherits(model, "biomod2_model") ){
    eval <- .evaluate.biomod2.formal.models(mod=model,data=data,stat=stat)
  } else {
    cat("\n\n! invalid model input => nothing returned !")
  }
  return(eval)
}

.evaluate.biomod2.formal.models <- function(mod, data, stat='TSS'){
  obs <- data[,mod@resp_name, drop=T]
  fit <- predict(mod, data[,mod@expl_var_names, drop=F], on_0_1000=T)
  if(stat != 'ROC'){
    thresh <- try(mod@model_evaluation[stat,'Cutoff'],silent=T)
  } else { thresh <- 500  } # no need to threshold
  if(inherits(thresh,"try-error")){
    thresh <- mod@model_evaluation[1,'Cutoff']
    cat("\n! no 'true' threshold defined for",stat,"... ", rownames(mod@model_evaluation)[1], "ones' taken.")
  }
  eval <- Find.Optim.Stat(Stat=stat, Fit=fit, Obs=obs, Fixed.thresh=thresh)
  return(eval)
}

.evaluate.biomod2.models.out <- function(mod, data, stat='TSS'){
  formal.models.names <- BIOMOD_LoadModels(mod)
  eval <- lapply(formal.models.names, function(x, data, stat){
    xx <- get(x)
    eval <- vapply(stat,
                   function(s){.evaluate.biomod2.formal.models(mod=xx, data=data, stat=s)},
                   FUN.VALUE = c(Evaluating.data=0, Cutoff=500, Sensitivity=0, Specificity=0))
    return(t(eval))
  }, data=data, stat=stat)
  names(eval) <- formal.models.names
  return(eval)
}

.evaluate.biomod2.ensemble.models.out <- function(mod, data, stat='TSS'){
  formal.models.names <- BIOMOD_LoadModels(mod)
  eval <- lapply(formal.models.names, function(x, data, stat){
    xx <- get(x)
    eval <- vapply(stat,
                   function(s){.evaluate.biomod2.formal.models(mod=xx, data=data, stat=s)},
                   FUN.VALUE = c(Evaluating.data=0, Cutoff=500, Sensitivity=0, Specificity=0))
    return(t(eval))
  }, data=data, stat=stat)
  names(eval) <- formal.models.names
  return(eval)
}

## TEST ##
# em <- evaluate(mod=mod.out,data=data,stat=c('TSS','ROC'))





######## LOW LEVEL FUCTIONS ##############

##' @name Find.Optim.Stat
##' @title Calculate the best score according to a given evaluation method
##'
##' @description \code{Find.Optim.Stat} is an internal \pkg{biomod2} function
##' to find the threshold to convert continuous values into binary ones leading
##' to the best score for a given evaluation metric.
##'
##' @usage
##'   Find.Optim.Stat(Stat='TSS',
##'                   Fit,
##'                   Obs,
##'                   Nb.thresh.test = 100,
##'                   Fixed.thresh = NULL)
##'
##' @param Stat either 'ROC', TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR',
##'         'POFD', 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR' or 'ORSS'
##' @param Fit vector of fitted values (continuous)
##' @param Obs vector of observed values (binary)
##' @param Nb.thresh.test integer, the numer of thresholds tested over the
##'        range of fitted value
##' @param Fixed.thresh integer, if not \code{NULL}, the only threshold value tested
##'
##' @details
##'   Please refer to \code{\link[biomod2]{BIOMOD_Modeling}} to get more information about this metrics.
##'   If you give a \code{Fixed.thresh}, no optimisation will be done. Only the score for this threshold will be returned.
##'
##' @return
##'   A 1 row x 4 column \code{matrix} :
##'   \itemize{
##'     \item{\code{best.iter}:}{ the best score obtained for chosen statistic}
##'     \item{\code{cutoff}:}{ the associated cut-off used for transform fitted vector into binary}
##'     \item{\code{sensibility}:}{ the sensibility with this threshold}
##'     \item{\code{specificity}:}{ the specificity with this threshold}
##'   }
##'
##' @author Damien Georges
##'
##' @seealso
##'   \code{\link[biomod2]{BIOMOD_Modeling}},
##'   \code{\link[biomod2]{getStatOptimValue}},
##'   \code{\link[biomod2]{calculate.stat}}
##'
##' @examples
##'   a <- sample(c(0,1),100, replace=TRUE)
##'
##'   ##' random drawing
##'   b <- runif(100,min=0,max=1000)
##'   Find.Optim.Stat(Stat='TSS',
##'                   Fit=b,
##'                   Obs=a)
##'
##'   ##' biased drawing
##'   BiasedDrawing <- function(x, m1=300, sd1=200, m2=700, sd2=200){
##'     return(ifelse(x<0.5, rnorm(1,m1,sd1), rnorm(1,m2,sd2)))
##'   }
##'
##'   c <- sapply(a,BiasedDrawing)
##'
##'   Find.Optim.Stat(Stat='TSS',
##'                   Fit=c,
##'                   Obs=a,
##'                   Nb.thresh.test = 100)
##'
##'
##' @keywords models
##' @keywords options
##' @keywords evaluation
##' 
Find.Optim.Stat <- 
  function(
    Stat = 'TSS',
    Fit,
    Obs,
    Nb.thresh.test = 100,
    Fixed.thresh = NULL
  ){

  ## remove all uninite values
  to_keep <- ( is.finite(Fit) & is.finite(Obs) )
  Fit <- Fit[to_keep]
  Obs <- Obs[to_keep]

  ## guess fit value scale (e.g. 0-1 for a classic fit or 0-1000 for a biomod2 model fit)
  fit.scale <- .guess.scale(Fit)

  ## check some data are still here.
  if(!length(Obs) | !length(Fit)){
    cat("Non finite obs or fit available => model evaluation skipped !")
    eval.out <- matrix(NA,1,4, dimnames = list(Stat, c("best.stat","cutoff","sensitivity","specificity")))
    return(eval.out)
  }

  if(length(unique(Obs)) == 1 | length(unique(Fit)) == 1){
#     warning("\nObserved or fited data contains only a value.. Evaluation Methods switched off\n",immediate.=T)
#     best.stat <- cutoff <- true.pos <- sensitivity <- true.neg <- specificity <- NA
      warning("\nObserved or fitted data contains a unique value... Be careful with this models predictions\n",immediate.=T)
      #best.stat <- cutoff <- true.pos <- sensitivity <- true.neg <- specificity <- NA
  } #else {
    if(Stat != 'ROC'){
      StatOptimum <- getStatOptimValue(Stat)
      if(is.null(Fixed.thresh)){ # test a range of threshold to get the one giving the best score
        if(length(unique(Fit)) == 1){
          valToTest <- unique(Fit)
          ## add 2 values to test based on maen with 0 and the guessed max of Fit (1 or 1000)
          valToTest <- round(c(mean(c(fit.scale["min"],valToTest)),
                               mean(c(fit.scale["max"],valToTest))))
        } else{
#           mini <- max(min(quantile(Fit,0.05, na.rm=T), na.rm=T),0)
#           maxi <- min(max(quantile(Fit,0.95, na.rm=T), na.rm=T),1000)
          mini <- max(min(Fit, na.rm=T), fit.scale["min"], na.rm = T)
          maxi <- min(max(Fit, na.rm=T), fit.scale["max"], na.rm = T)
          valToTest <- try(unique( round(c(seq(mini, maxi,
                                               length.out = Nb.thresh.test),
                                           mini, maxi)) ))
          if(inherits(valToTest, "try-error")){
            valToTest <- seq(fit.scale["min"], fit.scale["max"],
                             length.out = Nb.thresh.test)
          }
          # deal with unique value to test case
          if(length(valToTest)<3){
            valToTest <- round(c(mean(fit.scale["min"],mini), valToTest,
                                 mean(fit.scale["max"], maxi)))
          }
        }
#         valToTest <- unique( c(seq(mini,maxi,by=Precision), mini, maxi) )
      } else{
        valToTest <- Fixed.thresh
      }

      calcStat <- sapply(lapply(valToTest, function(x){return(table(Fit>x,Obs))} ), calculate.stat, stat=Stat)

      # scal on 0-1 ladder.. 1 is the best
      calcStat <- 1 - abs(StatOptimum - calcStat)

      best.stat <- max(calcStat, na.rm=T)

      cutoff <- median(valToTest[which(calcStat==best.stat)]) # if several values are selected

      misc <- table(Fit >= cutoff, Obs)
      misc <- .contagency.table.check(misc)
      true.pos <- misc['TRUE','1']
      true.neg <- misc['FALSE','0']
      specificity <- (true.neg * 100)/sum(misc[,'0'])
      sensitivity <- (true.pos * 100)/sum(misc[,'1'])
    } else{
      roc1 <- pROC::roc(Obs, Fit, percent=T, direction="<", levels = c(0,1))
      roc1.out <- pROC::coords(roc1, "best", ret = c("threshold", "sens", "spec"), transpose = TRUE)
      ## if two optimal values are returned keep only the first one
      if(!is.null(ncol(roc1.out))) roc1.out <- roc1.out[, 1]
      best.stat <- as.numeric(pROC::auc(roc1))/100
      cutoff <- as.numeric(roc1.out["threshold"])
      sensitivity <- as.numeric(roc1.out["sensitivity"])
      specificity <- as.numeric(roc1.out["specificity"])
    }
  #}
  eval.out <- cbind(best.stat,cutoff,sensitivity,specificity)
  rownames(eval.out) <- Stat

  return(eval.out)
}

getStatOptimValue <- function(stat){
  if(stat == 'TSS') return(1)
  if(stat == 'KAPPA') return(1)
  if(stat == 'PBC') return(1)
  if(stat == 'ACCURACY') return(1)
  if(stat == 'BIAS') return(1)
  if(stat == 'POD') return(1)
  if(stat == 'FAR') return(0)
  if(stat == 'POFD') return(0)
  if(stat == 'SR') return(1)
  if(stat == 'CSI') return(1)
  if(stat == 'ETS') return(1)
  if(stat == 'HK') return(1)
  if(stat == 'HSS') return(1)
  if(stat == 'OR') return(1000000)
  if(stat == 'ORSS') return(1)
}

calculate.stat <-
function(Misc, stat='TSS')
{
  # Contagency table checking
  Misc <- .contagency.table.check(Misc)

  # Defining Classification index
  hits <- Misc['TRUE','1']
  misses <- Misc['FALSE','1']
  false_alarms <- Misc['TRUE','0']
  correct_negatives <- Misc['FALSE','0']

  total <- sum(Misc)
  forecast_1 <- sum(Misc['TRUE',])
  forecast_0 <- sum(Misc['FALSE',])
  observed_1 <- sum(Misc[,'1'])
  observed_0 <- sum(Misc[,'0'])

  # Calculating choosen evaluating metric
  if(stat=='TSS'){
    return( (hits/(hits+misses)) + (correct_negatives/(false_alarms+correct_negatives)) -1 )
  }

  if(stat=='KAPPA'){
    Po <- (1/total) * (hits + correct_negatives)
    Pe <- ((1/total)^2) * ((forecast_1 * observed_1) + (forecast_0 * observed_0))
    return( (Po - Pe) / (1-Pe) )
  }

  if(stat=='ACCURACY'){
    return( (hits + correct_negatives) / total)
  }

  if(stat=='BIAS'){
    return( (hits + false_alarms) / (hits + misses))
  }

  if(stat=='POD'){
    return( hits / (hits + misses))
  }

  if(stat=='FAR'){
    return(false_alarms/(hits+false_alarms))
  }

  if(stat=='POFD'){
    return(false_alarms / (correct_negatives + false_alarms))
  }

  if(stat=='SR'){
    return(hits / (hits + false_alarms))
  }

  if(stat=='CSI'){
    return(hits/(hits+misses+false_alarms))
  }

  if(stat=='ETS'){
    hits_rand <- ((hits+misses)*(hits+false_alarms)) / total
    return( (hits-hits_rand) / (hits+misses+false_alarms-hits_rand))
  }

  if(stat=='HK'){
    return((hits/(hits+misses)) - (false_alarms/(false_alarms + correct_negatives)))
  }

  if(stat=='HSS'){
    expected_correct_rand <- (1/total) * ( ((hits+misses)*(hits+false_alarms)) +
      ((correct_negatives + misses)*(correct_negatives+false_alarms)) )
    return((hits+correct_negatives-expected_correct_rand) / (total - expected_correct_rand))
  }

  if(stat=='OR'){
    return((hits*correct_negatives)/(misses*false_alarms))
  }

  if(stat=='ORSS'){
    return((hits*correct_negatives - misses*false_alarms) / (hits*correct_negatives + misses*false_alarms))
  }

  if(stat=="BOYCE"){

  }

}

.contagency.table.check <- function(Misc){
  # Contagency table checking
  if(dim(Misc)[1]==1){
    if(row.names(Misc)[1]=="FALSE"){
      Misc <- rbind(Misc, c(0,0))
      rownames(Misc) <- c('FALSE','TRUE')
    } else{
      a <- Misc
    	Misc <- c(0,0)
  		Misc <- rbind(Misc, a)
      rownames(Misc) <- c('FALSE','TRUE')
  	}
  }

  if(ncol(Misc) != 2 | nrow(Misc) !=2 ){
    Misc = matrix(0, ncol=2, nrow=2, dimnames=list(c('FALSE','TRUE'), c('0','1')))
  }

  if((sum(colnames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) | (sum(rownames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) ){
    stop("Unavailable contingency table given")
  }

  if('0' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='0')] <- 'FALSE'
  if('1' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='1')] <- 'TRUE'

  return(Misc)
}


.guess.scale <- function(Fit){
  min <- 0
  max <- ifelse(max(Fit, na.rm = TRUE) <= 1, 1, 1000)
  out <- c(min, max)
  names(out) <- c("min", "max")
  return(out)
}
##' @include BiomodClass.R
##' @name FilteringTransformation
##' @aliases FilteringTransformation
##' @aliases FilteringTransformation-methods
##' @aliases FilteringTransformation,data.frame-method
##' @aliases FilteringTransformation,matrix-method
##' @aliases FilteringTransformation,numeric-method
##' @aliases FilteringTransformation,array-method
##' @aliases FilteringTransformation,RasterBrick-method
##' @aliases FilteringTransformation,RasterLayer-method
##' @aliases FilteringTransformation,RasterStack-method
##' 
##' @title Convert species' probability of occurrence into binary 
##' presence-absence data using a predefined threshold
##' 
##' @description
##' Function that converts an object containing probability values into 
##' a filtered object according to a pre-defined threshold(s).
##' 
##' 
##' @param data a numeric vector, a \code{matrix}, a \code{data.frame}, 
##' a \code{RasterLayer} or a \code{RasterStack} containing the data to 
##' be converted
##' @param threshold a numeric value or a vector containing the threshold
##' to be used for converting data.
##' 
##' @details
##' If data is a vector or a raster object, then the threshold should be a
##' numeric value. If data is matrix,dataframe or rasterStack, then the
##' threshold should have, in theory, as many values as the number of
##' columns or layers to transform.
##' In the particular case that the data to convert is a 
##' \code{matrix}/\code{data.frame} with several columns or a 
##' \code{RasterStack} with several layers and the threshold is a single
##' numeric value, the same threshold will be applied to all columns 
##' (resp. layers).  
##' 
##' @return 
##' An object of the same class than \code{data} with the values of data
##' if superior to \code{threshold} and 0 if not.
##' 
##' @author Wilfried Thuiller, Damien Georges
##'
##' @examples
##' xx <- rnorm(50,10)
##' yy <- FilteringTransformation(xx, 10)
##' 
##' cbind(xx,yy)
setGeneric("FilteringTransformation",
           function(data, threshold){
             standardGeneric("FilteringTransformation")
           })

setMethod('FilteringTransformation', signature(data='data.frame'),
  function(data, threshold)
  {
    data <- data.matrix(data)
    data[t(t(data)<threshold)] <-0

    ## check if some thresolds are NAs
    if(any(is.na(threshold))){
      data[,is.na(threshold)] <- NA
    }
    if(ncol(data)==1) data <- data[,1]
  	return(data)

  })

setMethod('FilteringTransformation', signature(data='matrix'),
  function(data, threshold)
  {
    data <- as.data.frame(data)
    return(FilteringTransformation(data, threshold))
  })

setMethod('FilteringTransformation', signature(data='numeric'),
  function(data, threshold)
  {
    data <- as.data.frame(data)
    return(FilteringTransformation(data, threshold))
  })

setMethod('FilteringTransformation', signature(data='array'),
          function(data, threshold)
          {
            if(length(dim(data)) == length(dim(threshold))){
              if(sum( dim(data)[-1] != dim(threshold)[-1] ) > 0 ){
                stop("data and threshold dimensions mismatch")
              }
            } else{
              if(sum( dim(data)[-1] != dim(threshold) ) > 0 ){
                stop("data and threshold dimensions mismatch")
              }
            }

            return(sweep(data,2:length(dim(data)),threshold,
                         function(x,y) {
                           if(!is.na(x)){
                             return(ifelse(x>y,x,0))
                           } else {
                             return(rep(NA,length(x)) )}
                         }))
          })


setMethod('FilteringTransformation', signature(data='RasterLayer'),
  function(data, threshold)
  {
    if(!is.na(threshold)){
      return(reclassify(data,c(-Inf,threshold,0)))
    } else{ ## return a empty map (NA everywhere)
      return(reclassify(data,c(-Inf,Inf,NA)))
    }
  })

setMethod('FilteringTransformation', signature(data='RasterStack'),
  function(data, threshold)
  {
    if(length(threshold) == 1){
      threshold <- rep(threshold, raster::nlayers(data))
    }
    StkTmp <- raster::stack()
    for(i in 1:raster::nlayers(data)){
      StkTmp <- raster::addLayer(StkTmp, FilteringTransformation(raster::subset(data,i,drop=TRUE), threshold[i]))
    }
    names(StkTmp) <- names(data)
    return(StkTmp)
  })

setMethod('FilteringTransformation', signature(data='RasterBrick'),
  function(data, threshold)
  {
    data <- raster::stack(data, RAT=FALSE)
    return(FilteringTransformation(data, threshold))
  })
`.functionkeep` <-
function(object, AIC)
{
    list(df.resid=object$df.resid, deviance=object$deviance, term=as.character(object$formula)[3], AIC=AIC)
}

`level.plot` <-
function(data.in, XY, color.gradient='red', cex=1, level.range=c(min(data.in),max(data.in)), show.scale=TRUE, title="level plot", SRC=FALSE, save.file="no", ImageSize="small", AddPresAbs=NULL, PresAbsSymbol=c(cex*0.8,16,4),...){

  extra.args <- list(...)
  if(!is.null(extra.args$multiple.plot)){
    multiple.plot <- extra.args$multiple.plot
  } else {
    multiple.plot <- FALSE
  }


    if(color.gradient!='grey' && color.gradient!='red' && color.gradient!='blue') stop("\n color.gradient should be one of 'grey', 'red' or 'blue' \n")
    if(ncol(XY)!=2) stop("\n wrong coordinates given in 'XY': there should be two columns \n")
    if(nrow(XY)!=length(data.in)) stop("\n data and coordinates should be of the same length \n")

#     if() multiple.plot <- TRUE  else multiple.plot <- FALSE


    if(SRC){
        if(length(unique(data.in))>4){
            cat("\n not possible to render SRC plot -> more than four different values in data ")
            SRC <- F
        } else {
            SRCvalues <- sort(unique(data.in))
			      col_id <- data.in + 3
            color.system <- c("red", "lightgreen", "grey", "darkgreen")
            title <- paste("SRC plot ", title, sep="")
        }
    } else
	{    if(color.gradient=='grey') {
#           color.system <- c()
#           for(i in seq(93,10,length.out=100)) color.system <- c(color.system, gray(i/100))
#           color.system <- c(gray(0.93), color.system, gray(0))
	        color.system <- gray(seq(0.95,0, length.out=102))

        }
        if(color.gradient=='blue') {
          color.system <- c('grey88',
          rainbow(45, start=0.5, end=0.65),
          rainbow(10, start=0.65, end=0.7),
          rainbow(45, start=0.7, end=0.85),
          'red')
        }
        if(color.gradient=='red') {
          color.system <- c(
          'grey88',
          c(rep(c(colors()[c(417,417,515)]), each=5),
          rev(rainbow(55, start=0.13, end=0.23 )),
          rev(rainbow(50, start=0.08, end=0.13 )[seq(1,50,length.out=15)]),
          rev(rainbow(50, end=0.08)[seq(1,50,length.out=15)])),
          'brown2')
        }


      #if range wanted is broader than possible, set to actual range limits
      #if(level.range[1]<min(data.in)) level.range[1] <- min(data.in)
      #if(level.range[2]>max(data.in)) level.range[2] <- max(data.in)

      #determine the color code to assess to each value

      # define a vector for make a correspundance between values and colors
      val_seq <- c(seq(level.range[1], level.range[2], length.out=101),Inf)
      col_id <- sapply(data.in, function(x){return(which(x<=val_seq)[1])})


#       g <- gg <- data.in
#       gg[gg <= level.range[1]]  <- level.range[1]
#       gg[gg >= level.range[2]] <- level.range[2]
#       gg <- gg-min(g)
#       gg <- gg/max(gg)*100 + 1
#
#       #over and under-ranged values set to limits of color range
#       gg[g < level.range[1]] <- 1
#       gg[g > level.range[2]] <- 102
    }

    # define plotting symbols for presence and absences if required by user
    if(!is.null(AddPresAbs)){
        cex2 <- PresAbsSymbol[1]
        pchPres <- PresAbsSymbol[2]
        pchAbs <- PresAbsSymbol[3]
    }

    #define image size for JPEG and TIFF
    if(ImageSize=="small") {SizeInPix <- 480; FontSize=12} else if(ImageSize=="standard") {SizeInPix <- 1000; FontSize=22} else if(ImageSize=="large") {SizeInPix <- 2000; FontSize=44}

    if(save.file == "pdf") pdf(paste(title, ".pdf", sep=""))
    if(save.file == "jpeg") jpeg(paste(title, ".jpeg", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize, quality=85)
    if(save.file == "tiff") tiff(paste(title, ".tiff", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize)
    if(save.file == "postscript") postscript(paste(title, ".eps", sep=""))

    if(show.scale){

        if(!multiple.plot) layout(matrix(c(1,2),nrow=1), widths=c(5,1), heights=c(1,1))
        plot(XY[,2]~XY[,1], col=color.system[col_id], cex=cex, pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)
        #Add Presence and Absence locations if requested by user:
        if(!is.null(AddPresAbs)){points(AddPresAbs[AddPresAbs[,3]==1,1:2], col="black", pch=pchPres, cex=cex2); points(AddPresAbs[AddPresAbs[,3]==0,1:2], col="black", pch=pchAbs, cex=cex2)}

        par(mar=c(0.1,0.1,0.1,0.1))
        plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
        polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)

        if(SRC){ legend(0,0.8,legend=list(' (1) new', ' (0) stable', '(-1) kept', '(-2) lost'),cex=1, fill=rev(color.system),bty='n')
         } else {
          if(level.range[1] == min(data.in)) lmin <- round(level.range[1], digits=2) else lmin <- paste(round(level.range[1], digits=2), " or lower", sep="")
          if(level.range[2] == max(data.in)) lmax <- round(level.range[2], digits=2) else {lmax <- paste(round(level.range[2], digits=2), " or over", sep="") }

          if(!multiple.plot){
              legend(0.2,0.92,legend=list(lmax,'','','','',round((3*level.range[2]+level.range[1])/4, digits=2),'','','','',round(sum(level.range)/2, digits=2),
              '','','','',round((level.range[2]+3*level.range[1])/4, digits=2),'','','','',lmin),cex=1, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')

          } else legend(0.2,1.05,legend=list(lmax,'','','','',round((3*level.range[2]+level.range[1])/4, digits=2),'','','','',round(sum(level.range)/2, digits=2),
          '','','','',round((level.range[2]+3*level.range[1])/4, digits=2),'','','','',lmin), cex=cex, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')

        }
    }
     else{
          plot(XY[,2]~XY[,1], col=color.system[col_id], cex=cex, pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)
          #Add Presence and Absence locations if requested by user:
          if(!is.null(AddPresAbs)){points(AddPresAbs[AddPresAbs[,3]==1,1:2], col="black", pch=pchPres, cex=cex2); points(AddPresAbs[AddPresAbs[,3]==0,1:2], col="black", pch=pchAbs, cex=cex2)}
     }

    if(!is.null(AddPresAbs)) rm(cex2,pchPres,pchAbs)
    if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()

}
##' @name makeFormula
##' @title Standardized formula maker
##' 
##' @description
##' makeFormula is an internal \pkg{biomod2} function that can be useful
##' to help users to build easily some standardized formula used later by
##' statistical models.
##' 
##' @param respName a \code{character} indicating the response variable
##' name
##' @param explVar a \code{matrix} or a \code{data.frame}, the 
##' explanatory variables table that will be considered at modelling step
##' @param type either 'simple', 'quadratic', 'polynomial' or 's_smoother'
##' defining the type of formula you want to build
##' @param interaction.level an \code{integer}, the interaction level
##' depth between explanatory variables
##' @param \ldots some additional arguments (see details)
##' 
##' @details
##' It is advised to give only a subset of \code{explVar} table to avoid
##' useless memory consuming. If some explanatory variables are factorial
##' ones, you have to give a \code{data.frame} for \code{explVar} where
##' associated columns are define as \code{factor}.
##' 
##' \code{...} argument available values are :
##' 
##' - `k` the smoothing parameter value (used only if 
##' \code{type = 's_smoother'}) corresponding to \code{k} parameter 
##' of \pkg{mgcv} \code{\link[mgcv]{s}}  or \code{df} \pkg{gam} 
##' \code{\link[gam]{s}} arguments.
##'
##' @return a \code{link[stats]{formula}} class object that can be
##' directly given to most of \R statistical models.
##' 
##' @author Damien Georges
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_ModelingOptions}}, 
##' \code{link[stats]{formula}}
##' 
##' @keywords models
##' @keywords formula
##' @keywords options
##' 
##' @examples
##' ##' create simulated data
##' myResp <- sample(c(0, 1), 20, replace = TRUE)
##' myExpl <- 
##'   matrix(
##'     runif(60), 
##'     ncol = 3, 
##'     dimnames=list(NULL, c('var1', 'var2', 'var3'))
##'   )
##' 
##' ##' create a formula
##' myFormula <- 
##'   makeFormula( 
##'     respName = 'myResp',
##'     explVar = head(myExpl),
##'     type = 'quadratic',
##'     interaction.level = 0
##'   )
##'   
##' ##' show formula created
##' myFormula
##'
makeFormula <- function(
  respName, 
  explVar, 
  type = 'simple', 
  interaction.level = 0, 
  ...
){
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  # This function return a string in a well formated way. May be give as formula argument to a "basic"
  # statistical model.
  # Several types of models are available
  #
  # D.GEORGES 12/2011
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  sup_args <- list(...)

  # 0. Supported Types =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  availableTypes = c("simple", "quadratic", "polynomial", "s_smoother")

  # 1. Check Given Args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(!is.character(respName) || length(respName)!=1){
    stop("You must give a unique response variable name")
  }

  if(!is.data.frame(explVar) &&  !is.matrix(explVar)){
    stop("You must give explanatory variable table")
  }

  if(!(type %in% availableTypes)){
    stop(paste("Formula type must be one of : ", toString(availableTypes), sep=""))
  }

  explVarNames <- colnames(explVar)
  if(respName %in% explVarNames){ # remove the response variable data if it's given
    explVar <- explVar[, - which(explVarNames == respName), drop=FALSE]
    explVarNames <- colnames(explVar)
  }

  interaction.level <- min(interaction.level, ncol(explVar))

  # 2. Create the formula =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  junk <- c(1)

  switch(EXPR=type,
         "simple" = {junk <- paste(junk,paste(explVarNames, collapse=" + "), sep=" + ") },

         "quadratic" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
#                junk <- c(junk, paste(explVarNames[v], "+I(", explVarNames[v],
#                                       "^2)+I(",explVarNames[v],"^3)", sep="") )
                junk <- paste(junk, paste(explVarNames[v], "+I(", explVarNames[v],
                                      "^2)", sep=""), sep=" + " )
             } else { junk <- paste(junk, explVarNames[v], sep=" + ") }
           } },

         "polynomial" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
               junk <- paste(junk, paste(explVarNames[v],
                                      "+I(", explVarNames[v],
                                      "^2)+I(",explVarNames[v],
                                      "^3)",sep=""), sep=" + " )
#                   junk <- c(junk, paste(explVarNames[v],
#                                       "+poly(",explVarNames[v],
#                                       ",3)",sep="") )
             } else { junk <- paste(junk, explVarNames[v], sep=" + ") }
           } },

         "s_smoother" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
               if(is.null(sup_args$k)){
                 junk <- paste(junk, paste("s(",explVarNames[v],
                                       ")",sep=""), sep=" + " )
               } else{
                 junk <- paste(junk, paste("s(",explVarNames[v],
                                       ",k=",sup_args$k,")",sep=""), sep=" + " )
               }

             } else { junk <- paste(junk, explVarNames[v], sep=" + ") }
           } })


  # interactions
  junk.inter <- NULL
  if(interaction.level > 0){
    for(i.l in 1:interaction.level){
      inter.tab <- combn(explVarNames,i.l+1)
      junk.inter <- paste(junk.inter, paste(apply(inter.tab,2,paste,collapse=":"),collapse=" + "), sep= " + ")
    }
  }

  if(length(junk.inter)){
    junk <- paste(junk,junk.inter,sep = "")
  }

  # 2. Return the formula =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  return(as.formula(paste(respName," ~ ", junk, sep="")))

}
############################
# models_scores_graph fct
# Damien G. - 2014/10/22
############################

## Description ##
# This function is a graphic tools to represent evaluation
# scores of models produced within biomod2 according to 2
# different evalution methods. Models can be grouped in several
# ways (by algo, by cv run, ...) to highlight potential differences
# in models quality due to chosen models, cross validation sampling 
# bias,...

## Input ##
# obj : an biomod2 modeling or ensemble-modeling object
# metrics : charcter vector of 2 chosen metrics (e.g c("ROC", "TSS"))
# by : the way models are grouped ('models', 'algos', 'cv_run' or 'data_set')
# plot : if you want to produce plot or not
# ... : several graphical options

## Ouput ##
# the ggplot2 object used to produce the graph is returned. That implies that 
# user should quite easily customize this plot.

## Main code ##
models_scores_graph <- function(obj, metrics = NULL, by = 'models', plot = TRUE, ...){
  
  ## get additional args
  args = list(...)
  
  ## check args
  ck_args <- .models_scores_graph_check_args(obj, metrics = metrics, by = by, args = args)
  
  metrics <- ck_args$metrics
  xlim <- ck_args$xlim
  ylim <- ck_args$ylim
  main <- ck_args$main
  
  ## get models scores
  scores <- get_evaluations(obj, as.data.frame = T )
  
  ## add some columns to enable different wy of grouping scores
  scores$mod = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[1] } )
  scores$alg = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[length(xx) - 2] } )
  scores$run = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[length(xx) - 1] } )
  scores$dat = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[length(xx)] } )
  
  ## extract some usefull info
  eval_data <- ifelse( all(is.na( scores$Evaluating.data )), "Testing.data", "Evaluating.data")
  
  ## calculate summaries statistics
  
  ### models mean scores calculation
  models_mean <- switch(by,
                        models = sapply( unique(scores$mod) , 
                                         function(x){ 
                                           sapply( metrics,  
                                                   function(xx, x){
                                                     mean(scores[ scores$mod == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                   },
                                                   x = x ) 
                                           }),
                        algos = sapply( unique(scores$alg) , 
                                         function(x){ 
                                           sapply( metrics,  
                                                   function(xx, x){
                                                     mean(scores[ scores$alg == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                   },
                                                   x = x ) 
                                         }),
                        cv_run = sapply( unique(scores$run) , 
                                      function(x){ 
                                        sapply( metrics,  
                                                function(xx, x){
                                                  mean(scores[ scores$run == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                },
                                                x = x ) 
                                      }),
                        data_set = sapply( unique(scores$dat) , 
                                           function(x){ 
                                             sapply( metrics,  
                                                     function(xx, x){
                                                       mean(scores[ scores$dat == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                     },
                                                     x = x ) 
                                           }) )
  
  ### sd of models scores calculation
  models_sd <- switch(by,
                        models = sapply( unique(scores$mod) , 
                                         function(x){ 
                                           sapply( metrics,  
                                                   function(xx, x){
                                                     sd(scores[ scores$mod == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                   },
                                                   x = x ) 
                                         }),
                        algos = sapply( unique(scores$alg) , 
                                       function(x){ 
                                         sapply( metrics,  
                                                 function(xx, x){
                                                   sd(scores[ scores$alg == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                 },
                                                 x = x ) 
                                       }),
                        cv_run = sapply( unique(scores$run) , 
                                      function(x){ 
                                        sapply( metrics,  
                                                function(xx, x){
                                                  sd(scores[ scores$run == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                },
                                                x = x ) 
                                      }),
                        data_set = sapply( unique(scores$dat) , 
                                           function(x){ 
                                             sapply( metrics,  
                                                     function(xx, x){
                                                       sd(scores[ scores$dat == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                     },
                                                     x = x ) 
                                           }) )
  
  ### merge data to fit with ggplot2 formalism
  ggdat <- merge( 
            data.frame(name = colnames(models_mean), t(models_mean) ),
            data.frame(name = colnames(models_sd), t(models_sd) ),
            by = "name" )
  
  colnames(ggdat) <- c("name", "mean1", "mean2", "sd1", "sd2")
  
  #   ggdat <- data.frame( 
  #             merge( 
  #               reshape2::melt(models_mean), 
  #               reshape2::melt(models_sd),
  #               by = c("X1","X2") ) )
  #   
  #   colnames(ggdat) <- c("metric","name","mean","sd")

  ## produce plots
  gg <- ggplot(ggdat, aes_string(x="mean1", y="mean2", colour="name", fill = NULL))
  
  ### add mean poins
  gg <- gg + geom_point()  

  
  ### add axis names and remove legend name
  gg <- gg + xlab(metrics[1]) + ylab(metrics[2]) + theme(legend.title=element_blank())
  
  ### fix scale
  if( length(ylim) | length(xlim)){
    gg <- gg + coord_cartesian(ylim=ylim, xlim=xlim)
  }
  
  ### add title
  if(length(main)){
    gg <- gg + labs(title=main)
  }
  
  ### add error bars
  limits1 <- aes_string(xmax = "mean1 + sd1", xmin= "mean1 - sd1", fill = NULL)
  limits2 <- aes_string(ymax = "mean2 + sd2", ymin= "mean2 - sd2", fill = NULL)
  gg <- gg + geom_errorbar(limits2,width=0) + geom_errorbarh(limits1, height=0)

  if(plot){
    print(gg)
  }
  
  invisible(gg)
} ## end of models_scores_graph function

.models_scores_graph_check_args <- function(obj, metrics = NULL, by = 'models', args){
  ## check obj type
  if(!inherits(obj, c("BIOMOD.models.out", "BIOMOD.EnsembleModeling.out"))){
    stop("obj must be a 'BIOMOD.models.out' or a 'BIOMOD.EnsembleModeling.out' object")
  }
  
  ## check metrics
  scores <- get_evaluations(obj, as.data.frame = T )
  
  avail_metrics <- unique( scores$Eval.metric )
  if( length(avail_metrics)<2 ){
    stop("at least 2 different evaluations metrics needed")
  }
  if (is.null(metrics)){
    metrics <- avail_metrics[1:2]
    warnings(toString(metrics), " evaluation metrics automatically selected")
  }
  
  ## check by
  if(! (by %in% c('models', 'algos', 'cv_run', 'data_set') ) ){
    stop("by arg should be one of 'models', 'algos', 'cv_run' or 'data_set'")
  }
  
  ## check extra args
  test_args_names <- ! ( names(args) %in% c("xlim", "ylim", "main", "col"))
  if( any ( test_args_names ) ){
    stop("unknown ", toString( names(args)[ test_args_names ] )," args")
  }
  
  
  xlim <- args$xlim
  ylim <- args$ylim
  main <- args$main
  
  return(list(metrics = metrics,
              xlim = xlim,
              ylim = ylim,
              main = main))
  
} ## end of checking args function

# x11()
# models_scores_graph(myBiomodModelOut, by = 'cv_run', metrics=c("ROC","TSS"))
# models_scores_graph(myBiomodEM, by = 'algos', metrics=c("ROC","TSS"))
# 
# str(gg)
# get_evaluations(myBiomodModelOut, as.data.frame=T)
`multiple.plot` <-
function(Data, coor, color.gradient='red', plots.per.window=9, cex=1, save.file="no", name="multiple plot", ImageSize="small", AddPresAbs=NULL, PresAbsSymbol=c(cex*0.8,16,4)){

    if(nrow(coor) != nrow(Data)) stop("Uncorrect mapping coordinates: coor and Data are not of the same length")
    if(color.gradient!='grey' && color.gradient!='red' && color.gradient!='blue') stop("\n color.gradient should be one of 'grey', 'red' or 'blue' \n")

#     assign("multiple", 564, pos=1) #to communicate to level.plot that a multiple plot is wanted (564 is just random)
    multiple <- 564

    #function plotting color boxes
    pbox <- function(co){
        plot(x=c(-1,1), y=c(0,1), xlim=c(0,1), ylim=c(0,1), type="n", axes=FALSE)
        polygon(x=c(-2,-2,2,2), y=c(-2,2,2,-2), col=co, border=NA)
    }

    #Take off NA data
    Data <- Data[ ,!is.na(Data[1,])]


    #calculating the number of windows to open
    NbPlots <- ncol(Data)
    NbWindows <- ceiling(NbPlots/plots.per.window)
    if(NbWindows==1) plots.per.window <- NbPlots

    #define image size for JPEG and TIFF
    if(ImageSize=="small") {SizeInPix <- 480; FontSize=12} else if(ImageSize=="standard") {SizeInPix <- 1000; FontSize=22} else if(ImageSize=="large") {SizeInPix <- 2000; FontSize=44}

    if(save.file=="pdf") pdf(paste(name, ".pdf", sep=""))
    if(save.file=="jpeg") jpeg(paste(name, ".jpeg", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize, quality=85)
    if(save.file=="tiff") tiff(paste(name, ".tiff", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize)
    if(save.file=="postscript") postscript(paste(name, ".eps", sep=""))

    for(W in 1:NbWindows){
        if(save.file=="no") dev.new()

        Wstart <- (W-1)*plots.per.window + 1
        if(W*plots.per.window > NbPlots) Wfinal <- NbPlots  else Wfinal <- W*plots.per.window
        DataW <- as.data.frame(Data[ ,Wstart:Wfinal])
        colnames(DataW) <- colnames(Data)[Wstart:Wfinal]

        #determine the organisation of the plots on the window
        W.width <- ceiling(sqrt(plots.per.window))
        W.height <- ceiling(plots.per.window/W.width)

        #create object for scaling the legend
        legendcex <- 0.64+1/exp(W.height)
#         assign("legendcex", 0.64+1/exp(W.height), pos=1)

        #matrix of indexes for ordering the layout
        mat <- c(1,2)
        for(i in 1:(W.width-1))  mat <- c(mat, mat[1:2] + 4*i)
        mat <- rbind(mat, mat+2)
        for(i in 1:(W.height-1))  mat <- rbind(mat, mat[1:2,] + W.width*4*i)

        layout(mat, widths=rep(c(1,0.3), W.width), heights=rep(c(0.2,1), W.height))

        par(mar = c(0.1,0.1,0.1,0.1))
        for(i in 1:(Wfinal-Wstart+1)){
             pbox("grey98")
             text(x=0.5, y=0.8, pos=1, cex=1.6, labels=colnames(DataW)[i], col="#4c57eb")
             pbox("grey98")
             level.plot(DataW[,i], XY=coor, color.gradient=color.gradient, cex=cex, title="", AddPresAbs=AddPresAbs, PresAbsSymbol=PresAbsSymbol, multiple.plot=TRUE)
        }

        #fill gaps by grey boxes
        if(W.width*W.height-plots.per.window != 0) for(i in 1:((W.width*W.height-plots.per.window)*4)) pbox("grey98")

    } #W loop

#     rm('legendcex', 'multiple', pos=1)
    if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()
}
.Prepare.Maxent.WorkDir <- function(Data, xy, calibLines=NULL, RunName=NULL,
                                    VarImport=0, evalData=NULL, evalxy=NULL,
                                    species.name=NULL, modeling.id='',
                                    background_data_dir = 'default'){
  cat('\n\tCreating Maxent Temp Proj Data..')

  ## initialise output
  MWD <- list()
  class(MWD) <- "maxent_workdir_info"

  ## default parameters setting
  if(is.null(RunName)) RunName <- colnames(Data)[1]
  if(is.null(species.name)) species.name <- colnames(Data)[1]
  if(is.null(calibLines)) calibLines <- rep(T,nrow(Data))

  ## define all paths to files needed by MAXENT.Phillips
  m_workdir <- file.path(species.name,'models',modeling.id,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
  # check wordir unicity
  while(file.exists(m_workdir)){
    m_workdir <- file.path(species.name,'models',modeling.id,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
  }
  MWD$m_workdir <- m_workdir

  m_outdir <- file.path(species.name,'models',modeling.id, paste(RunName,'_MAXENT.Phillips_outputs',sep=''))
  MWD$m_outdir <- m_outdir
  MWD$m_outputFile <- file.path(m_outdir,paste(RunName,'_Pred_swd.csv',sep=''))

  ## directories creation
  dir.create(m_workdir, showWarnings=FALSE, recursive=TRUE, mode='777')
  dir.create(m_outdir, showWarnings=FALSE, recursive=TRUE, mode='777')

  # Presences Data
  m_speciesFile <- file.path(m_workdir,"Sp_swd.csv")
  MWD$m_speciesFile <- m_speciesFile
  presLines <- which((Data[,1]==1) & calibLines)
  absLines <- which((Data[,1]==0) & calibLines)
  Sp_swd <- cbind(rep(RunName,length(presLines)),
                      xy[presLines,],
                      Data[presLines,2:ncol(Data),drop=FALSE])
  colnames(Sp_swd) <- c('species','X','Y',colnames(Data)[2:ncol(Data)])
  write.table(Sp_swd, file=m_speciesFile,  quote=FALSE, row.names=FALSE, sep=",")

  # Background Data
  ## create background file only if needed
  if(background_data_dir == 'default'){
    m_backgroundFile <- file.path(m_workdir,"Back_swd.csv")
    MWD$m_backgroundFile <- m_backgroundFile
    # keep only 0 of calib lines
    Back_swd <- cbind(rep("background",length(absLines)),xy[absLines,],Data[absLines,2:ncol(Data),drop=FALSE])
    colnames(Back_swd)  <- c("background",colnames(Back_swd)[-1])
    write.table(Back_swd, file=m_backgroundFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  } else { ## use background directory given as an option
    MWD$m_backgroundFile <- background_data_dir
  }

  # Prediction Data
  m_predictDir <- file.path(m_workdir,"Predictions")
  dir.create(m_predictDir, showWarnings=FALSE, recursive=TRUE)

  m_predictFile <- file.path(m_predictDir, "Pred_swd.csv")
  MWD$m_predictFile <- m_predictFile
  Pred_swd <- cbind(rep("predict",nrow(xy)),xy,Data[,2:ncol(Data),drop=FALSE])
  colnames(Pred_swd)  <- c("predict", colnames(xy), colnames(Data)[-1])
  write.table(Pred_swd, file=m_predictFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")


  ### USELESS ###
  # dealing with variable importances stuff
#   if( VarImport > 0){
#     for( vari in colnames(Data)[-1] )
#       for(vi in 1:VarImport){
#         proj_tmp <- Pred_swd
#         proj_tmp[,1] <- rep(paste(vari,'_',vi,sep=""),nrow(proj_tmp))
#         proj_tmp[,vari] <- sample(proj_tmp[,vari])
#         write.table(proj_tmp, file=file.path(species.name,'models',modeling.id,"MaxentTmpData","Pred",paste(vari,'_',vi,"_swd.csv",sep="")), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
#       }
#   }

  # dealing with independent evaluation data
  if(!is.null(evalData)){
    m_predictEvalFile <- file.path(m_predictDir, "PredEval_swd.csv")
    MWD$m_predictEvalFile <- m_predictEvalFile
    Pred_eval_swd <- cbind(rep("predictEval",nrow(evalxy)),evalxy,evalData[,2:ncol(evalData),drop=FALSE])
    colnames(Pred_eval_swd)  <- c("predict",colnames(Back_swd)[-1])
    write.table(Pred_eval_swd, file=m_predictEvalFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  }

  return(MWD)
}

.Delete.Maxent.WorkDir <- function(MWD, silent=FALSE){
  if(!silent) cat('\n\tRemoving Maxent Temp Data..')
  if(inherits(MWD, "maxent_workdir_info")){
    unlink(unique(sub("/part([[:digit:]]+)$", "", MWD$m_workdir)),recursive = TRUE, force = TRUE)
  } else{
    if(!silent) cat('\n\t! Invalid maxent work dir object -> MAXENT.Phillips temp files have not been removed')
  }
}

.create.maxent.bg.dir <- function(expl.stk, bm.dat, ...){
  args <-NULL
  args <- list(...)
  if(is.null(args$out.dir)) args$out.dir <- file.path(bm.dat@sp.name, "maxent.env.layers")
  if(is.null(args$NAflag)) args$NAflag <- -9999

  ## create the output directory
  dir.create(args$out.dir, showWarnings = FALSE, recursive = TRUE)
  ## save a copy of the environmental layers
  test.write <- sapply(names(bm.dat@data.env.var), function(ev){
    out.file <- file.path(args$out.dir, paste0(ev, ".asc"))
    cat("\n> writting", out.file)
    writeRaster(subset(expl.stk, ev),
                format = 'ascii',
                NAflag = args$NAflag,
                filename = out.file,
                overwrite = TRUE)
  })

  ## return the path to the sirectory where rasters are stored
  return(args$out.dir)
}

# Maxent Projection working directory preparation -=-=-=-=-=-=-=- #

setGeneric(".Prepare.Maxent.Proj.WorkDir",
            def = function(Data, ...){
              standardGeneric( ".Prepare.Maxent.Proj.WorkDir" )
            } )


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='data.frame'),
          def = function(Data, xy, species.name =".", proj.name=".", silent=FALSE){
            ## initialise output
            MWD <- list()
            class(MWD) <- "maxent_workdir_info"

            if(!silent) cat('\n\t\tCreating Maxent Temp Proj Data...')
            if(is.null(xy)) xy <- matrix(1,nrow=nrow(Data), ncol=2, dimnames=list(NULL, c("X","Y")))

            ## define all paths to files needed by MAXENT.Phillips
            m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            # check wordir unicity
            while(file.exists(m_workdir)){
              m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            }
            MWD$m_workdir <- m_workdir

            dir.create(m_workdir, recursive=TRUE, showWarnings=FALSE)

            # Proj Data
            m_predictFile <- file.path(m_workdir, "Pred_swd.csv")
            MWD$m_predictFile <- m_predictFile
            Proj_swd <- cbind(rep("proj",nrow(xy)),xy,Data)
            colnames(Proj_swd)  <- c("proj","X","Y",colnames(Data))
            write.table(Proj_swd, file=m_predictFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")

            return(MWD)
            })


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='RasterStack'),
          def = function(Data, species.name =".",proj.name=".", silent=FALSE, split.proj = 1){
            ## initialise output
            MWD <- list()
            class(MWD) <- "maxent_workdir_info"

            if(!silent) cat('\n\t\tCreating Maxent Temp Proj Data...')

            ## define all paths to files needed by MAXENT.Phillips
            m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            # check wordir unicity
            while(file.exists(m_workdir)){
              m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            }
#             MWD$m_workdir <- m_workdir
#             dir.create(m_workdir, recursive=TRUE, showWarnings=FALSE)

            ## create the list of extent our raster will be crop at
            pred.nrow <- nrow(Data)
            pred.ncol <- ncol(Data)
            seq.col <- round(seq(1, pred.ncol, length.out = split.proj + 1))
            ext.list <- lapply(1:split.proj, function(i){ return(extent(Data, 1, pred.nrow, seq.col[i], seq.col[i + 1]))})

            # Proj Data
            m_predictFile <- NULL
            for(spl in 1:split.proj){
              m_workdirTmp <- file.path(m_workdir, paste0("part", spl))
              dir.create(m_workdirTmp, showWarnings = FALSE, recursive = TRUE)
              MWD$m_workdir[[paste0("part", spl)]] <- m_workdirTmp
              for(l in names(Data)){
                m_predictFileTmp <- file.path(m_workdirTmp, paste0(l,'.asc'))
                if(! file.exists(m_predictFileTmp)){
                  if(!silent) cat("\n\t\t\t>",l ,"\t:\t" )
                  if(split.proj == 1){ ## no croping in this case => just write the raster as an ascii file
                    if(grepl(".asc", filename(raster::subset(Data,l,drop=TRUE)) ) ){
                      if(!silent) cat("copying ascii file")
                      file.copy(filename(raster::subset(Data,l,drop=TRUE)), m_predictFileTmp)
                    } else{
                      if(!silent) cat("creating ascii file")
                      writeRaster(raster::subset(Data,l,drop=TRUE), filename=m_predictFileTmp,
                                  format='ascii', overwrite=TRUE)
                    }
                  } else{ ## crop the raster within parts
                  crop(raster::subset(Data,l,drop=TRUE),
                       ext.list[[spl]], filename=m_predictFileTmp,
                       format='ascii', overwrite=TRUE)
                  }
                } else{
                  if(!silent) cat("\n", m_predictFileTmp ,'already created !')
                }
                m_predictFile <- c(m_predictFile, m_predictFileTmp)
              }
              MWD$m_predictFile[[paste0("part", spl)]] <- m_predictFile
            }
            return(MWD)
          })

setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='RasterLayer'),
          def = function(Data, species.name =".",proj.name=".", silent=FALSE){
            .Prepare.Maxent.Proj.WorkDir(Data = stack(Data), species.name = species.name,
                                         proj.name = proj.name , silent = silent)
          })

# .Prepare.Maxent.Proj.WorkDir <- function(Data, xy, proj_name=NULL){
#   cat('\n\tCreating Maxent Temp Proj Data..')
#
#   if(is.null(proj_name)) proj_name <- colnames(Data)[1]
#   dir.create(paste(getwd(),'/',proj_name,'/MaxentTmpData', sep=""), showWarnings=FALSE)
#
#   # Proj Data
#   Proj_swd <- cbind(rep("proj",nrow(xy)),xy,Data)
#   colnames(Proj_swd)  <- c("proj","X","Y",colnames(Data))
#   write.table(Proj_swd, file=paste(getwd(),'/',proj_name,"/MaxentTmpData/Proj_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
# }
#
# .Prepare.Maxent.Proj.Raster.WorkDir <- function(Data, proj_name=NULL){
#   cat('\n\tCreating Maxent Temp Proj Data..')
#
#   if(is.null(proj_name)){
#     stop("Please refere explicitly a proj name!")
#   }
#   dir.create(paste(getwd(),'/',proj_name,'/MaxentTmpData/Proj', sep=""), showWarnings=FALSE, recursive=TRUE)
#
#   # Proj Data
#   for(l in names(Data)){
#     if(! file.exists(file.path(proj_name,'MaxentTmpData','Proj',paste(l,'.asc',sep='')))){

#       if(grepl(".asc", filename(raster::subset(Data,l,drop=TRUE)) ) ){
#         cat("\n copying ascii file")
#         file.copy(filename(raster::subset(Data,l,drop=TRUE)), file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')))
#       } else{
#         cat("\n creating ascii file")
#         writeRaster(raster::subset(Data,l,drop=TRUE), filename=file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')),
#             format='ascii', overwrite=TRUE)
#       }
#
#     } else{
#       cat("\n", file.path(proj_name,'MaxentTmpData',paste(l,'.asc',sep='')),'already created !')
#     }
#
#   }
# }
##' @name ProbDensFunc
##' @title Probability Density Function
##' @description
##' Using a variety of parameters in modelling will inevitably bring
##' variability in predictions, especially when it comes to making future
##' predictions.
##' This function enables an overall viewing of the future predictions
##' range per species and gives the likelihood of range shift estimations.
##' It will calculate the optimal way for condensing a defined proportion
##' (50, 75, 90 and 95\% per default) of the data.
##' 
##' 
##' @param initial a vector (resp. a \code{SpatialPointsDataFrame}) in a
##' binary format (ones and zeros) representing the current distribution
##' of a species which will be used as a reference for the range change
##' calculations
##' @param projections a matrix (resp; a \code{rasterStack}) grouping all
##' the predictions where each column is a single prediction. Make sure
##' you keep projections in the same order as the initial vector (
##' line1 = site1, line2 = site2, etc.).
##' @param plothist set to TRUE to plot the range change histogram
##' @param cvsn stands for "current vs new". If true, the range change
##' calculations will be of two types: the percentage of cells currently
##' occupied by the species to be lost, and the relative percentage of
##' cells currently unoccupied but projected to be, namely 'new' cells,
##' compared to current surface range.
##' @param groups an option for ungrouping the projections enabling a
##' separated visualisation of the prediction range per given group. A
##' matrix is expected where each column is a single prediction and each
##' line is giving details of one parameter (See the examples section).
##' @param resolution the step used for classes of prediction in graphics.
##' The default value is 5.
##' @param filename the name of file (with extension) where plots will be
##' stored. If not \code{NULL}, no plotting windows will be open
##' @param ... further arguments:
##' 
##' - \code{lim}: ordered numeric vector indicating the proportion of data
##' to consider for histogram representation (by default : c(0.5, 0.75,
##' 0.9, 0.95))
##' - \code{nb.points.max}: the maximum number of points to sample, 25000
##' by default (useful for huge \code{raster*} objects)
##' 
##' @details
##' The future range changes are calculated as a percentage of the
##' species' present state. For example, if a species currently occupies
##' 100 cells and is estimated by a model to cover 120 cells in the 
##' future, the range change will be + 20\%.
##' 
##' Resolution : Note that modifying the resolution will directly
##' influence the probability scale. Bigger classes will accumulate a
##' greater number of predictions and therefore represent a greater
##' fraction of the total predictions. The probability is in fact that of
##' the class and not of isolated events.
##' 
##' @return This is a plotting function, no objects are returned or
##' created.
##' 
##' @author Wilfried Thuiller, Bruno Lafourcade, Damien Georges 
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' @keywords optimize
##' @keywords distribution
##' 
##' @examples 
##' \dontrun{
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"), row.names = 1)
##' head(DataSpecies)
##' 
##' ##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' ##' the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' ##' remove all 0 from response vector to work with
##' ##' presence only data (Pseudo Absences selections)
##' rm_id <- which(myResp==0)
##' myResp <- myResp[-rm_id]
##' 
##' 
##' ##' the XY coordinates of species data
##' myRespXY <- DataSpecies[-rm_id,c("X_WGS84","Y_WGS84")]
##' 
##' 
##' ##' Environmental variables extracted from BIOCLIM
##' myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio7.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio11.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio12.grd",
##'                              package="biomod2"))
##' 
##' ##' 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName,
##'                                      PA.nb.rep=3)
##' 
##' ##' 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' ##' 3. Doing Modelisation
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
##'                                      models = c('CTA','RF','GLM','GAM','ANN','MARS'),
##'                                      models.options = myBiomodOption,
##'                                      NbRunEval=5,
##'                                      DataSplit=70,
##'                                      Prevalence=0.5,
##'                                      models.eval.meth = c('TSS'),
##'                                      do.full.models = FALSE,
##'                                      rescal.all.models=T,
##'                                      modeling.id='test')
##' 
##' ##' 4. Build ensemble-models that will be taken as reference
##' myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
##'                                        chosen.models = 'all',
##'                                        em.by = 'all',
##'                                        eval.metric = c('TSS'),
##'                                        eval.metric.quality.threshold = c(0.7),
##'                                        prob.mean = TRUE,
##'                                        prob.median = TRUE)
##' 
##' ##' 5. Projection on future environmental conditions
##' 
##' ###' load future environmental conditions from biomod2 package
##' myExpl_fut <- raster::stack( system.file( "external/bioclim/future/bio3.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio4.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio7.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio11.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio12.grd",
##'                                   package="biomod2"))
##' 
##' myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
##'                                         new.env = myExpl_fut,
##'                                         proj.name = 'future',
##'                                         selected.models = 'all',
##'                                         binary.meth = 'TSS',
##'                                         compress = FALSE,
##'                                         build.clamping.mask = TRUE)
##' 
##' BIOMOD_EnsembleForecasting(projection.output=myBiomodProjection,
##'                            EM.output=myBiomodEM,
##'                            binary.meth='TSS')
##' 
##' ##' 6. load binary projections
##' consensusBin <- raster::stack('GuloGulo/proj_future/proj_future_GuloGulo_ensemble_TSSbin.grd')
##' projectionsBin <- raster::stack('GuloGulo/proj_future/proj_future_GuloGulo_TSSbin.grd')
##' 
##' ##' 7. build a ref state based on ensemble-models
##' ref <- sampleRandom(subset(consensusBin, 1, drop=T), size=5000, sp=T, na.rm=T)
##' 
##' ##' 8. autoatic creation of groups matrix
##' find_groups <- function(diff_by_pix){
##'   data.set <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='data.set')
##'   run.eval <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='run.eval')
##'   models <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='models')
##'   return(rbind(data.set,run.eval,models))
##' }
##' 
##' groups <- find_groups(projectionsBin)
##' 
##' ##' 9. plot ProbDensFunct graphs
##' ProbDensFunc(initial = ref,
##'              projections = projectionsBin,
##'              plothist=TRUE,
##'              cvsn=TRUE,
##'              groups=groups,
##'              resolution=2,
##'              filename=NULL,
##'              lim=c(0.5,0.8,0.95))
##' 
##' ###' 3 plots should be produced.. Should be convenient to save it within a device
##' ###' supporting multiple plots.
##' 
##' }
##' 
ProbDensFunc <- function(
  initial,
  projections,
  groups = NULL,
  plothist = TRUE,
  cvsn = TRUE,
  resolution = 5,
  filename = NULL,
  ...
){

  args <- .ProbDensFunc.checkArgs(
    initial,
    projections,
    groups,
    plothist,
    cvsn,
    resolution,
    filename,
    ...
  )

  initial = args$initial
  projections = args$projections
  groups = args$groups
  plothist = args$plothist
  cvsn = args$cvsn
  resolution = args$resolution
  filename = args$filename
  lim = args$lim
  one_plot = args$one_plot

  rm(list='args')

  # results will be stored in out
  out <- list()

  # oppen the divice if filename
  if(length(filename)){
    switch(tools::file_ext(filename),
           pdf = pdf(filename),
           jpeg = jpeg(filename),
           tiff = tiff(filename),
           eps = postscript(filename),
           png = png(filename))
  }

  # area stores the species range change calculations
  area <- (apply(projections,2,sum, na.rm=T) / sum(initial, na.rm=T) - 1 ) * 100
  a <- round( (min(area, na.rm=TRUE)-(resolution+10))/10 ) *10
  b <- round( (max(area, na.rm=TRUE)+(resolution+10))/10 ) *10
  area_hist <- hist(area, breaks = seq(a,b,resolution), plot=FALSE)
  area_hist$density <- area_hist$counts / sum(area_hist$counts)


  # analysis of the distribution density and calculation of the probability of events
  area_sorted <- sort(area)
  nb <- round(length(area_sorted) * lim)
  lower_limit <- upper_limit <- c()

  for(i in 1:length(lim)){
     g <- rep(NA,length(area_sorted)-nb[i])
     for(j in 1:length(g)) g[j] <- diff(range(area_sorted[j:min(j+nb[i],length(area_sorted))]))
     lower_limit <- c(lower_limit,area_sorted[which.min(g)])
     upper_limit <- c(upper_limit,area_sorted[which.min(g)]+g[which.min(g)])
  }

  names(lower_limit) <- names(upper_limit) <- paste(lim*100, "%", sep="")

  out$stats <- cbind(lower_limit,upper_limit)

  if(!is.null(groups)){
    if(! (length(filename) | one_plot)) dev.new()
    par(mfrow=c(1,nrow(groups)))

    color.samp <- list()
    for(pa in 1:nrow(groups)){

      lv <- levels(as.factor(as.matrix(groups[pa,])))
      color.samp[[pa]] <- colors()[sample(c(90,417,552,616,382,11,150,468,28,31,420,476,333),length(lv))]
      g <- hist(area, breaks = seq(a,b,resolution), plot=FALSE)
      fac <- (max(g$counts) / sum(g$counts)) / max(g$density)
      g$density <- g$density * fac

      plot(g, freq=FALSE, border='grey88', main=row.names(groups)[pa], xlab="Species range change (%)", ylab="Event   occurence   probability")
      lines(density(area, width=30)$x, density(area, width=30)$y*fac)
      for(i in 1:length(lv)){
        div <- length(area) / length(area[groups[pa,]==lv[i]])
        lines(density(area[groups[pa,]==lv[i]],width=30)$x, density(area[groups[pa,]==lv[i]], width=30)$y /div*fac, col=color.samp[[pa]][i])
      }
      lv <- as.factor(as.matrix(groups[pa,]))
      leg <- list()
      for(j in 1:length(levels(lv))) leg[[j]] <- levels(lv)[j]
      legend("topright", legend=leg, bty='n',fill=color.samp[[pa]])
    }
  }

  if(cvsn){
    #calculation of the 2 axes independently (lost vs new sites)
    area2 <- (apply(projections[which(initial==1),],2,sum) / sum(initial==1) -1) * 100
    area3 <- area - area2

    if(! (length(filename) | one_plot)) dev.new()
    par(mfrow=c(1,nrow(groups)))

    for(i in 1:nrow(groups)){
      lv <- as.factor(as.matrix(groups[i,]))
      leg <- list()
      for(j in 1:length(levels(lv))) leg[[j]] <- levels(lv)[j]

      levels(lv) <- 1:length(levels(lv))
      plot(area3~area2, xlab='current', ylab='new', ylim=c(0,if(max(area3)<100){100}else{max(area3)+30}), xlim=c(-100,0), col=color.samp[[i]][lv], pch=20, main=row.names(groups)[i])
      legend("bottomleft", legend=leg, bty='n',fill=color.samp[[i]])
      abline(0,-1, col='grey80')
      abline(100,-1, col='grey80')
      text(x=-97,y=103,pos=1,label="SRC = 0", col="black", cex=0.8)
      if(max(area3)<100)text(x=-3,y=103,pos=1,label="SRC = 100", col="black", cex=0.8)
      else text(x=if(max(area3)<200){-(max(area3)+33-100)}else{-96},y=if(max(area3)<200){max(area3)+33}else{203},pos=1,label="SRC = 100", col="black", cex=0.8)
    }
  }

  #if plot of distribution plot wanted
  if(plothist){
    if(! (length(filename) | one_plot)) dev.new()
    par(mfrow=c(1,1))

    col_list <- colorRampPalette(c("dodgerblue1","steelblue1","slategray1","aliceblue"))(length(lim))

    for( l in length(lim):1){
      hist( mean(out$stat[l,], na.rm=T), breaks=out$stat[l,], col=col_list[l], xlim=c(a,b), ylim=c(0,max(area_hist$density)*1.2), xlab="",ylab="", main="", add=ifelse(l==length(lim),FALSE,TRUE))
      abline(v=out$stat[l,],col= col_list[l], lwd=1.7)
    }

    legend("topright",
          legend=rownames(out$stat), bty='n', fill=col_list, cex=0.8, title='distrib. of data' )

    par(new=TRUE)
    plot(area_hist, freq=FALSE, col="white", xlim=c(a,b), ylim=c(0,max(area_hist$density)*1.2), main="Probability density function", xlab="Species range change (%)", ylab="Event occurence probability")
  }

  if(length(filename)) dev.off()

  return(out)
}

.ProbDensFunc.checkArgs <- function(initial,
                                    projections,
                                    groups,
                                    plothist,
                                    cvsn,
                                    resolution,
                                    filename,
                                    ...){
  add.args <- list(...)
  if(is.null(add.args$nb.points.max) | !is.numeric(add.args$nb.points.max)) add.args$nb.points.max <- 25000
  if(is.null(add.args$lim) | !is.numeric(add.args$lim)) add.args$lim <- c(0.5, 0.75, 0.90, 0.95)
  if(is.null(add.args$one_plot)) add.args$one_plot <- FALSE

  # check lim arg
  if(sum(add.args$lim>1 | add.args$lim<0) > 0 ) stop("'lim' must be a numeric vector with 0 to 1 values")


  # check args types
  if(inherits(projections, 'Raster')){
    if(!inherits(initial, c('RasterLayer', 'SpatialPointsDataFrame')))
      stop("If projections is a raster object, initial should be a 'RasterLayer' or a 'SaptialPointDataFrame'")
  } else if(is.matrix(projections)){
    if(!is.numeric(initial)){
      stop("If projections is a matrix, initial should be a 'numeric'")
    }
  } else{
    stop("projections should be a 'matrix' or a 'RasterStack'")
  }

  # extract values
  if(inherits(projections, 'Raster')){
    if(inherits(initial, 'SpatialPointsDataFrame')){
      if(nrow(initial) > add.args$nb.points.max){
        initial[sort(sample(1:nrow(initial),size=add.args$nb.points.max)), drop=FALSE]
      }
    } else {
      initial <- sampleRandom(initial, size=min(add.args$nb.points.max, ncell(initial)) ,sp=TRUE, na.rm=TRUE)
    }

    projections <- extract(projections, initial, method='simple', na.rm=FALSE)
    initial <- initial@data[,1]
  } else{
    if(length(initial) != nrow(projections))
      stop("initial & projections dimensions don't match")
    if(length(initial) > add.args$nb.points.max){
      kept_rows <- sort(sample(1:length(initial), size=add.args$nb.points.max, replace=FALSE))
      initial <- initial[kept_rows]
      projections <- projections[kept_rows,,drop=FALSE]
    }
  }

  # remove NAs
  na_rows <- unique(c( which(is.na(initial)), which(is.na(projections), arr.ind=TRUE)[,1] ) )
  if(length(na_rows)){
    initial <- initial[-na_rows]
    projections <- projections[-na_rows,,drop=FALSE]
  }

  # check groups arg
  if(is.null(groups)){
    if(cvsn){
      cat("\n\t! 'cvsn' was automatically switch off because no 'groups' given")
       cvsn <- FALSE
    }
  } else {
    if(!is.matrix(groups)) stop("'groups' should be a matrix")
    if(ncol(groups)!=ncol(projections)) stop("'groups' and 'projections' do not have the same number of columns (resp. layers)")
  }

  # check saving options
  if(!is.null(filename)){
    if( ! (tools::file_ext(filename) %in% c("pdf","jpeg","tiff","eps","png"))){
      filename <- paste(tools::file_path_sans_ext(filename),".pdf",sep="")
      cat("\n\t! 'filename' extension unknown => outputs will be stored in :", filename)
    }
  }

  return(list( initial = initial,
               projections = projections,
               groups = groups,
               plothist = plothist,
               cvsn = cvsn,
               resolution = resolution,
               filename = filename,
               lim = add.args$lim,
               one_plot = add.args$one_plot))
}

setGeneric( "Projection",
            def = function(models.name,
                           modeling.work.dir = getwd(),
                           new.env.data ,
                           ...){
                            standardGeneric( "Projection" )
                            } )

setMethod( 'Projection', signature(new.env.data = 'data.frame'),
  function(models.name,
           modeling.work.dir = getwd(),
           new.env.data ,
           xy = NULL,
           proj.name = NULL,
           binary.proj = NULL,
           filtred.proj = NULL,
           models.evaluation = NULL,
           models.options = NULL,
           compress = TRUE,
           scaled.models=TRUE,
           do.stack = FALSE){

#     # 1. loading resuired libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#     .Models.dependencies(silent=TRUE, models.options=models.options )

    # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

    sp.name  <- .extractModelNamesInfo(model.names=models.name, info='species')
    PA.run   <- .extractModelNamesInfo(model.names=models.name, info='data.set')
    eval.run <- .extractModelNamesInfo(model.names=models.name, info='run.eval')
    algo.run <- .extractModelNamesInfo(model.names=models.name, info='models')

    # 3. Printing Projection Summary =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

    # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#     cat('\nDoing Models Projections...')
#     if(length(grep('EF.',models.name)) > 0 ){
#       kept.models.name <- models.name[-grep('EF.',models.name)]
#       kept.algo.run <- algo.run[-grep('EF.',algo.run)]
#     } else {
#       kept.models.name <- models.name
#       kept.algo.run <- algo.run
#     }

    proj.array <- lapply(models.name, .Projection.do.proj, env=new.env.data, xy=xy, scaled.models=scaled.models, proj.name=paste("proj_",proj.name, sep=""), models.options=models.options)
    proj.array <- as.data.frame(proj.array)
    names(proj.array) <- models.name

    # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(binary.proj)>0){
      cat("\nBinary transformations...")
      lapply(binary.proj, function(bin.proj){

        cuts <- unlist(lapply(names(proj.array), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
          }))

        proj.bin.array <- BinaryTransformation(proj.array, cuts)
        proj.bin.array <- DF_to_ARRAY(proj.bin.array)

        eval(parse(text = paste(proj.name,"_",sp.name,"_bin_",bin.proj, "_array <- proj.bin.array", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_bin_",bin.proj,
                                "_array, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_bin_",bin.proj,"_array' )",sep="")))

        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_bin_",bin.proj,"_array , proj.bin.array, cuts)", sep="" )))
      })
    }


    # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(filtred.proj)>0){
      cat("\nFiltered transformations...")
      lapply(filtred.proj, function(filt.proj){

        cuts <- unlist(lapply(names(proj.array), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
        }))

        proj.filt.array <- FilteringTransformation(proj.array, cuts)
        proj.filt.array <- DF_to_ARRAY(proj.filt.array)

        eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_array <- proj.filt.array", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                "_array, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_filt_",filt.proj,"_array' )",sep="")))

        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_array , proj.filt.array, cuts)", sep="" )))
      })
    }

    proj.array <- DF_to_ARRAY(proj.array)

    # 7. Saving projection on hard disk
    eval(parse(text = paste(proj.name,"_",sp.name, " <- proj.array", sep="")))
    eval(parse(text = paste("save(",proj.name,"_",sp.name, ", file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                            proj.name,"_",sp.name,"' )",sep="")))
    eval(parse(text = paste("rm(",proj.name,"_",sp.name,")", sep="" )))
    gc(reset=TRUE)

    return(invisible(proj.array))
  })








setMethod( 'Projection', signature(new.env.data = 'RasterStack'),
  function(models.name,
           modeling.work.dir = getwd(),
           new.env.data ,
           xy = NULL,
           proj.name = NULL,
           binary.proj = NULL,
           filtred.proj = NULL,
           models.evaluation = NULL,
           models.options = NULL,
           stack = TRUE,
           compress = TRUE,
           scaled.models=TRUE,
           do.stack = FALSE){

    # 1. loading resuired libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
    .Models.dependencies(silent=TRUE, models.options=models.options)

    # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

    sp.name  <- .extractModelNamesInfo(model.names=models.name, info='species')
    PA.run   <- .extractModelNamesInfo(model.names=models.name, info='data.set')
    eval.run <- .extractModelNamesInfo(model.names=models.name, info='run.eval')
    algo.run <- .extractModelNamesInfo(model.names=models.name, info='models')


    # 3. Printing Projection Summary =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

    # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    cat('\nDoing Models Projections...')
#     if(length(grep('EF.',models.name)) > 0 ){
#       kept.models.name <- models.name[-grep('EF.',models.name)]
#       kept.algo.run <- algo.run[-grep('EF.',algo.run)]
#     } else {
#       kept.models.name <- models.name
#       kept.algo.run <- algo.run
#     }

    if(do.stack){
      proj.ras <- lapply(models.name, .Projection.do.proj, env=new.env.data, scaled.models=scaled.models, proj.name=paste("proj_",proj.name, sep=""), models.options=models.options)

      # transform list of rasterLayers into a rasterStack
      proj.stack <- stack(x = proj.ras)

      names(proj.stack) <- models.name #names(proj.ras.mod)
      rm(proj.ras)

      # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(length(binary.proj)>0){
        cat("\nBinary transformations...")
        lapply(binary.proj, function(bin.proj){

          cuts <- unlist(lapply(names(proj.stack), function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
            }))

          proj.bin.stack <- BinaryTransformation(proj.stack, cuts)
          names(proj.bin.stack) <- paste(names(proj.stack), ".bin", sep="")

          eval(parse(text = paste(proj.name,"_",sp.name,"_bin_",bin.proj, "_RasterStack <- proj.bin.stack", sep="")))
          eval(parse(text = paste("save(",proj.name,"_",sp.name,"_bin_",bin.proj,
                                  "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                  proj.name,"_",sp.name,"_bin_",bin.proj,"_RasterStack' )",sep="")))

          eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_bin_",bin.proj,"_RasterStack , proj.bin.stack, cuts)", sep="" )))
        })
      }

      # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(length(filtred.proj)>0){
        cat("\nFiltered transformations...")
        lapply(filtred.proj, function(filt.proj){

          cuts <- unlist(lapply(names(proj.stack), function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
          }))

          proj.filt.stack <- FilteringTransformation(proj.stack, cuts)
          names(proj.filt.stack) <- paste(names(proj.stack), ".filt", sep="")

          eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_RasterStack <- proj.filt.stack", sep="")))
          eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                  "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                  proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterStack' )",sep="")))

          eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterStack , proj.filt.stack, cuts)", sep="" )))
        })
      }

      # 7. Saving projection on hard disk
       eval(parse(text = paste(proj.name,"_",sp.name, "_RasterStack <- proj.stack", sep="")))
       eval(parse(text = paste("save(",proj.name,"_",sp.name, "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                               proj.name,"_",sp.name,"_RasterStack' )",sep="")))
       eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_RasterStack)", sep="" )))
       gc(reset=TRUE)
    } else{

      # all models will be saved separatly
      proj.stack <- c() # list of saved files
      for(m.n in models.name){

        proj.ras <- .Projection.do.proj(m.n, env=new.env.data, scaled.models=scaled.models, proj.name=paste("proj_",proj.name, sep=""), models.options=models.options)
        names(proj.ras) <- m.n #names(proj.ras.mod)

        # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(length(binary.proj)>0){
          cat("\nBinary transformations...")
          lapply(binary.proj, function(bin.proj){

            cuts <- unlist(lapply(names(proj.ras), function(x){
              mod <- tail(unlist(strsplit(x,"_")), 3)[3]
              run <- tail(unlist(strsplit(x,"_")), 3)[2]
              dat <- tail(unlist(strsplit(x,"_")), 3)[1]
              return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
              }))

            proj.bin.ras <- BinaryTransformation(proj.ras, cuts)
            names(proj.bin.ras) <- paste(names(proj.ras), ".bin", sep="")

            eval(parse(text = paste(proj.name,"_",m.n,"_bin_",bin.proj, "_RasterLayer <- proj.bin.ras", sep="")))
            eval(parse(text = paste("save(",proj.name,"_",m.n,"_bin_",bin.proj,
                                    "_RasterLayer, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                    proj.name,"_",m.n,"_bin_",bin.proj,"_RasterLayer' )",sep="")))

            eval(parse(text = paste("rm(",proj.name,"_",m.n,"_bin_",bin.proj,"_RasterLayer , proj.bin.ras, cuts)", sep="" )))
          })
        }

        # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(length(filtred.proj)>0){
          cat("\nFiltered transformations...")
          lapply(filtred.proj, function(filt.proj){

            cuts <- unlist(lapply(names(proj.ras), function(x){
              mod <- tail(unlist(strsplit(x,"_")), 3)[3]
              run <- tail(unlist(strsplit(x,"_")), 3)[2]
              dat <- tail(unlist(strsplit(x,"_")), 3)[1]
              return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
            }))

            proj.filt.ras <- FilteringTransformation(proj.ras, cuts)
            names(proj.filt.ras) <- paste(names(proj.ras), ".filt", sep="")

            eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_RasterLayer <- proj.filt.ras", sep="")))
            eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                    "_RasterLayer, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                    proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterLayer' )",sep="")))

            eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterLayer , proj.filt.ras, cuts)", sep="" )))
          })
        }

        # 7. Saving projection on hard disk
        eval(parse(text = paste(proj.name,"_",m.n, "_RasterLayer <- proj.ras", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",m.n, "_RasterLayer, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                 proj.name,"_",m.n,"_RasterLayer' )",sep="")))
        proj.stack <- c( proj.stack, paste(proj.name,"_",m.n,"_RasterLayer",sep="") )
        eval(parse(text = paste("rm(",proj.name,"_",m.n,"_RasterLayer)", sep="" )))
        gc(reset=TRUE)
      }

    }

    ## remove MAXENT.Phillips tmp dir if exists
    if(file.exists(file.path(sp.name, proj.name, 'MaxentTmpData'))){
      .Delete.Maxent.WorkDir( file.path(sp.name, proj.name) )
    }

    return(invisible(proj.stack))
  })


setGeneric( ".Projection.do.proj",
            def = function(model.name, env, model.dir = NULL,...){
                    standardGeneric( ".Projection.do.proj" )
                    } )

setMethod('.Projection.do.proj', signature(env='data.frame'),
  function(model.name, env, xy = NULL, model.dir = NULL, scaled.models=TRUE, proj.name=NULL, models.options=NULL){
    cat('\n\t>', model.name)
    # automaticly fill model.dir if not given
    if(is.null(model.dir)){
      model.dir <- paste(getwd(),'/',.extractModelNamesInfo(model.name, info='species'),'/models', sep="")
    }

    # loading model
    if(length(c(grep('SRE',model.name) )) == 0){
      model.sp = eval(parse(text = load(paste(model.dir,'/',model.name, sep=""))) )
      eval(parse(text=paste("rm(",model.name,")",sep="")))
    }

    # check model.type
    model.type <- tail(unlist(strsplit(model.name, split="_")),1)
    if(!( model.type %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF', 'MAXENT.Phillips') )){
      if(!grep('EF.',model.type))
        stop('Unknown model type')
    }

    if(model.type == 'ANN'){
      set.seed(555) # to be able to refind our trees MAY BE BAD
      # proj automaticly scaled
      return(data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env, type = "raw")),
                                              name = model.name ) * 1000)))
    }

    if(model.type == 'CTA'){
      proj <- as.integer(as.numeric(predict(model.sp, env,type="prob")[,2]) * 1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }

      return( data.frame( proj = proj) )
    }

    if(model.type == 'FDA'){
      # proj automaticly scaled
      return( data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env, type = "posterior")[, 2]),
                                              name = model.name) * 1000)))
    }

    if(model.type == 'GBM'){
      best.iter <- gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
      proj <- as.integer(predict.gbm(model.sp, env, best.iter, type = "response") * 1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }

      return(data.frame( proj = proj))
    }

    if(model.type == 'GLM'){
      proj <- as.integer(.testnull(model.sp, Prev=0.5, env) * 1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }

      return( data.frame(proj = proj) )
    }

    if(model.type == 'GAM'){
      proj <- as.integer(.testnull(model.sp, Prev=0.5, env) * 1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }

      return( data.frame( proj = proj ) )
    }

#     if(model.type == 'MARS'){
#       # proj automaticly scaled
#       return(data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env)),
#                                               name = model.name) * 1000)))
#     }

    if(model.type == 'RF'){
      proj <- as.integer(as.numeric(predict(model.sp,env, type='prob')[,'1']) *1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }

      return( data.frame( proj = proj ))
    }

    if(model.type == 'SRE'){
      # loading data of the correspunding run
      load(paste(model.dir,'/',model.name,'/Data_',model.name, sep=""))
      return(eval(parse(text=paste("sre(Data_",model.name,"$Response, Data_",
                                   model.name,"$Explanatory, env, Data_",model.name,
                                   "$Quant)*1000", sep=""))))
    }

    if(model.type == 'MAXENT.Phillips'){
      if(!is.null(xy)){
        proj <- as.integer(predict( object=model.sp, newdata=env, proj_name=proj.name, xy=xy) * 1000)

#         if(scaled.models){
          proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
#         }

        return( data.frame( proj = proj ))
#
#         .Prepare.Maxent.Proj.WorkDir(env, xy, proj.name=file.path(.extractModelNamesInfo(model.name, info='species'), proj.name ))
#
#         cat("\t Running Maxent...")
#
#         system(command=paste("java -cp ", file.path(models.options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"), " density.Project \"", model.dir,.Platform$file.sep,
#                              model.name, .Platform$file.sep ,sub("_MAXENT.Phillips","",model.name),
#                              ".lambdas\" ", .extractModelNamesInfo(model.name, info='species'), "/", proj.name, "/MaxentTmpData/Proj_swd.csv ", .extractModelNamesInfo(model.name, info='species'), "/", proj.name, "/MaxentTmpData/projMaxent", sep=""), wait = TRUE)
#
#         maxent.proj <- read.asciigrid(paste(species.name=.extractModelNamesInfo(model.name, info='species'), "/", proj.name , "/MaxentTmpData/projMaxent.asc", sep=""))@data
#         .Delete.Maxent.WorkDir(species.name=paste(species.name=.extractModelNamesInfo(model.name, info='species'), "/", proj.name,sep=""))
#         return(proj = as.integer(.Rescaler5(as.numeric(maxent.proj[,1]),
#                                               name = model.name) * 1000))
      } else {
        cat('\n MAXENT.Phillips need coordinates to run! NA returned ')
        return(data.frame(rep(NA,nrow(env))))
      }
    }

  })








setMethod('.Projection.do.proj', signature(env='RasterStack'),
  function(model.name, env, model.dir = NULL, scaled.models=TRUE, proj.name=NULL, models.options=NULL){
    cat('\n\t>', model.name)

    # automaticly fill model.dir if not given
    if(is.null(model.dir)){
      model.dir <- paste(getwd(),'/',.extractModelNamesInfo(model.name, info='species'),'/models', sep="")
    }

    # loading model
    if(length(c(grep('SRE',model.name))) == 0){
      model.sp = eval(parse(text = load(paste(model.dir,'/',model.name, sep=""))) )
      eval(parse(text=paste("rm(",model.name,")",sep="")))
    }

    # check model.type
    model.type <- tail(unlist(strsplit(model.name, split="_")),1)
    if(!( model.type %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF', 'MAXENT.Phillips') )){
      if(!grep('EF.',model.type))
        stop('Unknown model type')
    }

    if(model.type == 'ANN'){
      set.seed(555) # to be able to refind our trees MAY BE BAD
      proj.ras <- predict(env, model.sp, type="raw")
      proj.ras[!is.na(proj.ras[])] <- .Rescaler5(proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(proj.ras * 1000))
    }

    if(model.type == 'CTA'){
      set.seed(123) # to be able to refind our trees MAY BE BAD
      proj.ras <- predict(env, model=model.sp, type='prob', index=2)
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )
    }

    if(model.type == 'FDA'){
      pred.ras <- predict(env, model.sp, type="post", index=2)
      pred.ras[!is.na(pred.ras[])] <- .Rescaler5(pred.ras[!is.na(pred.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(pred.ras * 1000))
    }

    if(model.type == 'GBM'){
      if(file.exists(paste(model.dir,'/',model.name,'_best.iter'))){
        load(paste(model.dir,'/',model.name,'_best.iter'))
      } else{
        best.iter <- gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
      }

      proj.ras <- predict(env, model.sp, n.trees=best.iter, type='response')
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }

      return( round(proj.ras*1000) )
    }

    if(model.type == 'GLM'){
      proj.ras <- predict(env, model=model.sp, type='response')
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }

      return( round(proj.ras*1000) )
    }

    if(model.type == 'GAM'){
      proj.ras <- predict(env, model=model.sp, type='response')
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )
    }

    if(model.type == 'MARS'){
      pred.ras <- predict(env, model.sp)
      pred.ras[!is.na(pred.ras[])] <- .Rescaler5(pred.ras[!is.na(pred.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(pred.ras * 1000) )
    }

    if(model.type == 'RF'){
      proj.ras <- predict(env, model=model.sp, type='prob', index=2)
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )
    }

    if(model.type == 'SRE'){
#       cat('\n SRE prediction not supported yet ! ')
      load(paste(model.dir,'/',model.name,'/Data_',model.name, sep=""))
      data.sre <- get(paste('Data_',model.name, sep=""))
      rm(list=paste('Data_',model.name, sep=""))
#       sre.out <- eval(parse(text=paste("sre(Data_",model.name,"$Response, Data_",
#                                    model.name,"$Explanatory, env, Data_",model.name,
#                                    "$Quant)*1000", sep="")))
      sre.out <- raster::subset(sre(data.sre$Response, data.sre$Explanatory, env, data.sre$Quant), 1, drop=TRUE) * 1000

      return(sre.out)
    }

    if(model.type == 'MAXENT.Phillips'){
      proj.ras <- predict( object=model.sp, newdata=env, proj_name=proj.name)

#       if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5(proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                   name=model.name, original=FALSE)
#       }

      return(round(proj.ras*1000))



    }
  })

.pseudo.absences.sampling <-
function(sp, env, nb.repet=1, strategy='random', distMin=0, distMax=NULL, nb.points=NULL, quant.SRE = 0, PA.table = NULL){

  # 1. Parameters checking
  args <- .check.params.pseudo.absences.sampling(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE)

  sp <- args$sp
  env <- args$env
  nb.repet <- args$nb.repet
  strategy <- args$strategy
  distMin <- args$distMin
  distMax <- args$distMax
  nb.points <- args$nb.points
  quant.SRE <- args$quant.SRE

  rm("args")

  if( (nb.repet == 0 | nb.points <= 0) & strategy != 'user.defined'){
    out <- NULL
  } else {
    out <- switch(strategy,
                   user.defined = user.defined.pseudo.abs.selection(sp, env, PA.table),
                   random = random.pseudo.abs.selection( sp, env, nb.points, nb.repet ),
                   sre = sre.pseudo.abs.selection(sp, env, quant.SRE, nb.points, nb.repet),
                   disk = disk.pseudo.abs.selection(sp, env, distMin, distMax, nb.points, nb.repet))
  }

  return(out)

#   # 2. Check if NA are present in sp or not to determine which dataset to use
#   if(sum(is.na(sp@data)) > 0 ){ # PA will be taken into response variable
#     cat("\n*** PA selection")
#     pa.tab <- switch(strategy,
#                      random = random.pseudo.abs.selection(data=sp, nb.points=nb.points, nb.repet=nb.repet),
#                      sre = sre.pseudo.abs.selection(sp),
#                      disk = disk.pseudo.abs.selection(sp))
#     .arranging.pa.table()
#   } else{ # PA will be taken into explanatory variables
#     if(inherits(env, 'Raster')){ # Raster env var case
#
#     } else if(inherits(env, 'SpatialPoints')){ # spatial data.frame case
#
#     }
#
#   }
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.check.params.pseudo.absences.sampling <- function(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE){
  cat("\n   > Pseudo Absences Selection checkings...")

  # define here the implemented strategies
  availableStrategies <- c("random", "sre", "disk", "user.defined")

  # 1. sp input checking
  if(is.vector(sp)){
    sp <- SpatialPointsDataFrame(matrix(0,ncol=2,nrow=length(sp)), data.frame(sp))
  }

  if(!(inherits(sp, 'SpatialPoints'))){
    stop("species input must be a SpatialPointsDataFrame object")
  }

  # 2. env input checking
  if(is.matrix(env) | is.data.frame(env)){
    if(nrow(env) != nrow(sp)){
      stop("Species and Explanatory must have same dimensions")
    }
    env <- SpatialPointsDataFrame(coordinates(sp), as.data.frame(env))
  }

  if(!inherits(env, 'SpatialPoints') & !inherits(env, 'Raster')){
    stop("Explanatory variables input must be a SpatialPointsDataFrame or a RasterStack object")
  }

  # 3. Strategy checking
  if( ! (strategy %in% c("random", "sre", "user.defined")) ){
    if( ( sum(abs(coordinates(sp))) == 0 ) | !( strategy %in% availableStrategies ) ){ # no coordinates or unknow strategy
      strategy <- "random"
      cat("\n   ! Random strategy was automatically selected (that can be due to points coordinates lack or unavailable strategy choosen)")
    }
  }

  # 4. Nb points checking
  if(strategy != "user.defined"){
    if(is.null(nb.points)){
      stop("You must give the number of pseudo absences you want")
    } else{
      nbTrueAbs <- .get.nb.true.abs(sp)
      if(nbTrueAbs >= nb.points){
        cat("\n    ! There is more 'true absences' than desired pseudo absences. No pseudo absences selection done.")
        nb.points = 0
        #       #### Return a flag that tell to function that no PA selected
        #       return(NULL)
      } else {
        nb.points = nb.points - nbTrueAbs
      }
    }
  }


  # 4. Nb repetition checking

  # 5. Distances checking
  if(!is.null(distMin)){
    if(distMin < 0){
        distMin <- 0
    }
  }

  if(!is.null(distMax)){
    if(distMax < 0){
        distMax <- NULL
    }
  }

  if(!is.null(distMax) & !is.null(distMin)){
    if(distMin >= distMax){
      stop("distMin >= distMax")
    }
  }

  # 6. SRE quantil checking
  if(strategy == 'SRE'){
    if( quant.SRE >= 0.5 | quant.SRE <0 ){
      stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
    }
  }

  # 7. return checked params
  return(list(sp = sp,
              env = env,
              nb.repet = nb.repet,
              strategy = strategy,
              distMin = distMin,
              distMax = distMax,
              nb.points = nb.points,
              quant.SRE = quant.SRE))

}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.get.nb.true.abs <- function(sp){
  if(is.vector(sp)) return(sum(sp==0, na.rm=TRUE))

  if(inherits(sp, 'SpatialPoints')) return(sum(sp@data==0, na.rm=TRUE))

  if(inherits(sp, 'Raster')) return(sum(sp[]==0, na.rm=TRUE))
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.nb.available.pa.cells <- function(data, PA.flag = NA){
  if(is.vector(data)){
    return(ifelse(is.na(PA.flag), sum(is.na(data)), sum(data == PA.flag, na.rm = TRUE)))
  }
  if(is.data.frame(data) | is.matrix(data)){
    return(ifelse(is.na(PA.flag), sum(is.na(data)), sum(data == PA.flag, na.rm = TRUE)))
  }
  if(inherits(data, 'SpatialPoints')){
    return(ifelse(is.na(PA.flag), sum(is.na(data@data)), sum(data@data == PA.flag, na.rm = TRUE)))
  }
  if(inherits(data, 'Raster')){
    return(ifelse(is.na(PA.flag), sum(is.na(data[])), sum(data[] == PA.flag, na.rm = TRUE)))
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.rand.pseudo.abs.selection <- function(data, nb.points){
  if(is.vector(data)){ return(sample(1:length(data), nb.points, replace=FALSE)) }

  if(inherits(data, 'SpatialPoints')){ return(sample(1:nrow(data@data), nb.points, replace=FALSE))}

  if(inherits(data, 'Raster')){ return(sort(sampleRandom(x=data, size=nb.points, cells=T)[,"cell"]))}
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
  setGeneric( "random.pseudo.abs.selection",
              def = function(sp,env, ...){
                      standardGeneric( "random.pseudo.abs.selection" )
                      } )
# }

setMethod('random.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function( sp, env, nb.points, nb.repet ){
            cat("\n   > random pseudo absences selection")

            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(sp)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(is.na(sp@data))
              for(j in 1:ncol(pa.tab)){
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(x = as.data.frame(env),
                                                         mask.out = pa.tab[, j, drop = FALSE])
                if(length(fact.level.cells)){
                  pa.tab[fact.level.cells, j] <- TRUE
                  cand.cells <- setdiff(cand.cells, fact.level.cells)
                }
                pa.tab[sample(x = cand.cells,
                              size = nb.points - length(fact.level.cells),
                              replace = FALSE), j] <- TRUE
              }
              return(list(xy = coordinates(sp),
                          sp = as.vector(sp@data),
                          env = as.data.frame(env@data),
                          pa.tab = pa.tab))
            } else {
              cat("\nUnsupported case yet!")
              return(NULL)
            }
          })

setMethod('random.pseudo.abs.selection', signature(env="RasterStack"),
          function( sp, env, nb.points, nb.repet ){
#             require('raster',quietly=T)
            cat("\n   > random pseudo absences selection")

            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(sp)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(is.na(sp@data))
              for(j in 1:ncol(pa.tab)){
                pa.tab[sample(x=cand.cells,size=nb.points,replace=FALSE),j] <- TRUE
              }
              env <- as.data.frame(extract(env, coordinates(sp)))

              return(list(xy = coordinates(sp),
                          sp = as.numeric(unlist(sp@data, use.names=FALSE)),
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              # create a mask containing all not already sampled points (presences and absences)
              mask.env <- mask.out <- raster::subset(env, 1, drop = TRUE)
              mask.env <- raster::reclassify(mask.env, c(-Inf, Inf, -1)) ## the area we want to sample
              mask.out[] <- NA

              # add presences and true absences in our mask
              in.cells <- cellFromXY(mask.env, coordinates(sp))
              mask.env[in.cells] <- NA
              mask.out[in.cells] <- 1

              # checking of nb candidates
              nb.cells <- .nb.available.pa.cells(mask.env, PA.flag = -1)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }

              # select cells into raster
              pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
              for( j in 1:ncol(pa.tab.tmp)){
                ## initialise the vector of sample cells
                SR <- NULL
                ## define a compy of the sampling mask
                mask.env.tmp <- mask.env
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(env, mask.out = mask.out)
                if(length(fact.level.cells)){
                  SR <- c(SR, fact.level.cells)
                  ## update the mask by removing already selected cells
                  mask.env.tmp[SR] <- NA
                }
                SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                         size = nb.points - length(SR),
                                         cells = T,
                                         na.rm = T)[, "cell", drop = TRUE])
                ## repeat sampling until haing the right number of points
                while(length(SR)<nb.points){
                  ## update the mask by removing already selected cells
                  mask.env.tmp[SR] <- NA
                  ## extract the missing number of points
                  SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                           size = nb.points - length(SR),
                                           cells = T,
                                           na.rm = T)[, "cell", drop = TRUE])
                }
                pa.tab.tmp[,j] <- SR
              }

              # puting cells in good format
              selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              for( j in 1:ncol(pa.tab)){
                pa.tab[is.element(selected.cells, pa.tab.tmp[,j]), j] <- TRUE
              }

              # puting presences, true absences and pseudo absences together
              xy <- rbind(coordinates(sp), xyFromCell(mask.env, selected.cells))
              xy <- .add_PA_rownames(xy)
              sp <- as.numeric(unlist(c(as.vector(sp@data),
                                        rep(NA, length(selected.cells))),
                                      use.names = FALSE))
              env <- extract(env, xy)

              pa.tab <- rbind(matrix(TRUE, nrow = (nrow(xy) - length(selected.cells)),
                                     ncol = ncol(pa.tab)), pa.tab)

              return(list(xy = xy,
                          sp = sp,
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))

            }
          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
setGeneric( "user.defined.pseudo.abs.selection",
            def = function(sp,env, ...){
              standardGeneric( "user.defined.pseudo.abs.selection" )
            } )
# }

setMethod('user.defined.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function( sp, env, pa.table ){
            cat("\n   > User defined pseudo absences selection")

              return(list(xy = coordinates(sp),
                          sp = as.vector(sp@data),
                          env = as.data.frame(env@data),
                          pa.tab = pa.table))

          })

setMethod('user.defined.pseudo.abs.selection', signature(env="RasterStack"),
          function( sp, env, pa.table ){
#             require('raster',quietly=T)
            cat("\n   > User defined pseudo absences selection")

#             env <- as.data.frame(extract(env, coordinates(sp), method='bilinear'))
              env <- as.data.frame(extract(env, coordinates(sp)))


            return(list(xy = coordinates(sp),
                        sp = as.numeric(unlist(sp@data, use.names=FALSE)),
                        env = as.data.frame(env),
                        pa.tab = pa.table))

          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
  setGeneric( "sre.pseudo.abs.selection",
              def = function(sp,env, ...){
                      standardGeneric( "sre.pseudo.abs.selection" )
                      } )
# }
setMethod('sre.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function(sp, env, quant.SRE, nb.points, nb.repet){
            cat("\n   > SRE pseudo absences selection")

            # 1. calculate sre to determine availables
            mask.in <- sre(Response = sp, Explanatory = env, NewData = env@data, Quant = quant.SRE)
            ## we want to sample PA out of the SRE => have to revert the mask
            mask.in <- data.frame(mask.in = !as.logical(mask.in))
#             # mask of already sampled points (presneces/absences)
#             mask[mask == 0] <- NA
#             mask[which(as.vector(sp@data)==1),1] <- 1
#             mask[which(as.vector(sp@data)==0),1] <- 0

            # 2. Check if NA are present in sp or not to determine which dataset to use
#             if(.nb.available.pa.cells(mask) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(mask.in$mask.in, PA.flag = TRUE)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(!mask.in$mask.in)
              for(j in 1:ncol(pa.tab)){
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(as.data.frame(env),
                                                         mask.out = pa.tab[, j, drop = FALSE],
                                                         mask.in = mask.in)
                pa.tab[c(fact.level.cells,
                         sample(x = setdiff(cand.cells, fact.level.cells),
                                size = nb.points - length(fact.level.cells),
                                replace = FALSE)), j] <- TRUE
              }
              return(list(xy = coordinates(sp),
                          sp = as.vector(sp@data),
                          env = as.data.frame(env@data),
                          pa.tab = pa.tab))

#             }
          })



setMethod('sre.pseudo.abs.selection', signature(env="RasterStack"),
          function(sp, env, quant.SRE, nb.points, nb.repet){
            cat("\n   > SRE pseudo absences selection")

            # 1. calculate sre to determine availables
            ## mask in which we want to sample
            mask.in <- sre(Response = sp, Explanatory = env, NewData = env, Quant = quant.SRE)
            ## remove all points that are in the cpecies SRE
            mask.in[mask.in[] > 0] <- NA

            ## mask of already sampled points (presences/absences)
            mask.out <- raster::subset(env, 1)
            mask.out[] <- NA; mask.out[cellFromXY(mask.out,
                                                  coordinates(sp)[is.element(as.vector(sp@data), c(0, 1)), ])]

#             # removing cells in envelops, presences and absences
#             mask[mask == 1] <- NA
#             mask[cellFromXY(mask, coordinates(sp)[which(as.vector(sp@data) ==1 ), ])] <- NA
#             mask[cellFromXY(mask, coordinates(sp)[which(as.vector(sp@data) == 0), ])] <- NA


            # checking of nb candidates
            nb.cells <- .nb.available.pa.cells(mask.in, PA.flag = 0)
            if(nb.cells <= nb.points){
              nb.repet <- 1
              nb.points <- nb.cells
              cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
            }

            # select cells into raster
            pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
            for( j in 1:ncol(pa.tab.tmp)){
              ## initialise the vector of sample cells
              SR <- NULL
              ## define a compy of the sampling mask
              mask.in.tmp <- mask.in
              ## force to get at least one value of each factorial variable
              fact.level.cells <- sample.factor.levels(env,
                                                       mask.out  = mask.out,
                                                       mask.in = mask.in)
              if(length(fact.level.cells)){
                SR <- c(SR, fact.level.cells)
                ## update the mask by removing already selected cells
                mask.in.tmp[SR] <- NA
              }
              SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                       size = nb.points - length(SR),
                                       cells = TRUE,
                                       na.rm = TRUE)[, "cell", drop = TRUE])
              ## repeat sampling until haing the right number of points
              while(length(SR) < nb.points){
                ## update the mask by removing already selected cells
                mask.env.tmp[SR] <- NA
                ## extract the missing number of points
                SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                         size = nb.points - length(SR),
                                         cells = T,
                                         na.rm = T)[, "cell", drop = TRUE])
              }
              pa.tab.tmp[, j] <- SR
            }

            # puting cells in good format
            selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
            pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
            colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
            for( j in 1:ncol(pa.tab)){
              pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
            }

            # puting presences, true absences and pseudo absences together
            xy <- rbind(coordinates(sp)[which(!is.na(as.vector(sp@data))),],
                        xyFromCell(mask.in, selected.cells))
            xy <- .add_PA_rownames(xy)
            sp <- as.numeric(unlist(c(na.omit(as.vector(sp@data)), rep(NA,length(selected.cells))), use.names=FALSE))
            env <- extract(env, xy)

            pa.tab <- rbind(matrix(TRUE,nrow=(nrow(xy)-length(selected.cells)), ncol=ncol(pa.tab)),
                           pa.tab)

            return(list(xy = xy,
                        sp = sp,
                        env = as.data.frame(env),
                        pa.tab = as.data.frame(pa.tab)))

          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setGeneric( "disk.pseudo.abs.selection",
              def = function(sp,env, ...){
                      standardGeneric( "disk.pseudo.abs.selection" )
                      } )

setMethod('disk.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function(sp, env, distMin, distMax, nb.points, nb.repet){
            cat("\n   > Disk pseudo absences selection")

            # 1. determining selectable area
            coor <- coordinates(sp)
            pres <- which(sp@data[,1]==1)
            true.abs <- which(sp@data[,1]==0)
            tmp.abs <- which(is.na(sp@data[,1])) #(1:ncol(sp@data))[-c(pres,true.abs)]
            outside <- rep(0, length(abs))
            inside <- rep(0, length(abs))


            for(i in 1:length(pres)){
              # removing points too close from presences
              inside <- inside + ( sqrt((coor[tmp.abs,1]-coor[pres[i],1])^2 + (coor[tmp.abs,2]-coor[pres[i],2])^2) > distMin )
              # keeping points not to far from presences
              if(!is.null(distMax)){
                outside <- outside + ( sqrt((coor[tmp.abs,1]-coor[pres[i],1])^2 + (coor[tmp.abs,2]-coor[pres[i],2])^2) < distMax )
              }
            }
            if(is.null(distMax)){ # no cells are too far
              outside <- outside + 1
            }
            selected.abs <- tmp.abs[ (inside == length(pres)) & (outside > 0) ]

            # 2. adding presences and true absences and selecting randomly pseudo absences

            return(random.pseudo.abs.selection( sp[c(pres, true.abs, selected.abs),],
                                                env[c(pres, true.abs, selected.abs),],
                                                nb.points, nb.repet ))


          })


## TODO(damien): remimplement disk.pseudo.abs.selection to call random.pseudo.abs.selection
setMethod('disk.pseudo.abs.selection', signature(env="RasterStack"),
          function(sp, env, distMin, distMax, nb.points, nb.repet){
            cat("\n   > Disk pseudo absences selection")

            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              env.tmp <- SpatialPointsDataFrame(coords = coordinates(sp),
#                                                 data = as.data.frame(extract(env,coordinates(sp),method='bilinear'))
                                                data = as.data.frame(extract(env,coordinates(sp))))

              return(disk.pseudo.abs.selection(sp, env.tmp, distMin, distMax, nb.points, nb.repet))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")

              # create a mask
              dist.mask <- raster::subset(env,1, drop=TRUE)
              dist.mask[] <- NA

              pres.xy <- coordinates(sp[which(sp@data[,1]==1),])
              dist.mask[cellFromXY(dist.mask,pres.xy)] <- 1

              dist.mask <- raster::distance(dist.mask)
              dist.mask <- mask(dist.mask, raster::subset(env,1, drop=TRUE))

              if(is.null(distMax)) distMax <- Inf

              mask.in <- reclassify(dist.mask, c(-Inf,distMin,NA ,distMin, distMax,-1, distMax,Inf,NA))

              # get the mask of already sampled mask
              mask.out <- mask.in; mask.out[] <- NA; mask.out[cellFromXY(mask.out, coordinates(sp))] <- 1

              # checking of nb candidates
              nb.cells <- .nb.available.pa.cells(mask.in, PA.flag = -1)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }

              # select cells into raster
              pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
              for( j in 1:ncol(pa.tab.tmp)){
                ## initialise the vector of sample cells
                SR <- NULL
                ## define a compy of the sampling mask
                mask.in.tmp <- mask.in
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(env, mask.out = mask.out)
                if(length(fact.level.cells)){
                  SR <- c(SR, fact.level.cells)
                  ## update the mask by removing already selected cells
                  mask.in.tmp[SR] <- NA
                }
                SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                         size = nb.points - length(SR),
                                         cells = T,
                                         na.rm = T)[, "cell", drop = TRUE])
                ## repeat sampling until haing the right number of points
                while(length(SR)<nb.points){
                  ## update the mask by removing already selected cells
                  mask.in.tmp[SR] <- NA
                  ## extract the missing number of points
                  SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                           size = nb.points - length(SR),
                                           cells = T,
                                           na.rm = T)[, "cell", drop = TRUE])
                }
                pa.tab.tmp[,j] <- SR
              }

              # puting cells in good format
              selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              for( j in 1:ncol(pa.tab)){
                pa.tab[is.element(selected.cells, pa.tab.tmp[,j]), j] <- TRUE
              }

              # puting presences, true absences and pseudo absences together
              xy <- rbind(coordinates(sp), xyFromCell(mask.in,
                                                      selected.cells))
              xy <- .add_PA_rownames(xy)
              sp <- as.numeric(unlist(c(as.vector(sp@data), rep(NA,length(selected.cells))),
                                      use.names = FALSE))
              env <- extract(env, xy)

              pa.tab <- rbind(matrix(TRUE,
                                     nrow = (nrow(xy) - length(selected.cells)),
                                     ncol = ncol(pa.tab)), pa.tab)

              return(list(xy = xy,
                          sp = sp,
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))

            }
          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# automaticaly add rownames to a data.frame
.add_PA_rownames <- function(xy){
  rn <- row.names(xy)
  missing_rn <- which(rn == "")
  if(length(missing_rn)){
    rn[missing_rn] <- paste("pa", 1:length(missing_rn), sep="")
  }
  rownames(xy) <- rn
  return(xy)
}
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# setMethod('disk.pseudo.abs.selection', signature(env="RasterStack"),
#           function(sp, env, distMin, distMax, nb.points, nb.repet){
#             cat("\n   > Disk pseudo absences selection")
#
#             # 1. Check if NA are present in sp or not to determine which dataset to use
#             if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
#               env.tmp <- SpatialPointsDataFrame(coords = coordinates(sp),
#                                                 data = as.data.frame(extract(env,coordinates(sp),method='bilinear')))
#
#               return(disk.pseudo.abs.selection(sp, env.tmp, distMin, distMax, nb.points, nb.repet))
#             } else {
#               cat("\n   > Pseudo absences are selected in explanatory variables")
#
#               # create a mask
#               mask <- maskInside <- maskOutside <- reclassify(raster::subset(env,1), c(-Inf,Inf,0))
#               pres.xy <- coordinates(sp[which(sp@data[,1]==1),])
#
#               # to convert longitudinal degrees into metters
#               coef.conversion <- ifelse(grepl("longlat",env@crs@projargs), 111319.5, 1)
#               #               coef.conversion <- 1
#               ## progress bar
#               cat("\n")
#               pb <- txtProgressBar(min = 0, max = nrow(pres.xy), initial = 0, char = "=-",width = 20,  style = 3, file = "")
#               for(i in 1:nrow(pres.xy)){
#                 setTxtProgressBar(pb,i)
#                 if(distMin > 0){
#                   maskInside <- maskInside + (distanceFromPoints(mask, pres.xy[i,]) > (distMin * coef.conversion))
#                 }
#                 if(!is.null(distMax)){
#                   maskOutside <- maskOutside + (distanceFromPoints(mask, pres.xy[i,]) <= (distMax * coef.conversion))
#                 }
#               }
#
#               if(distMin > 0){
#                 maskInside <- maskInside == nrow(pres.xy)
#               } else { # keep all cells
#                 maskInside <- maskInside + 1
#               }
#
#               if(!is.null(distMax)){
#                 maskOutside <- maskOutside > 0
#               } else{ # keep all cells
#                 maskOutside <- maskOutside + 1
#               }
#
#
#               mask <- maskInside * maskOutside
#               mask[mask==0] <- NA
#               mask <- (-1) * mask
#
#               # remove presences and true absences from our raster
#               mask[cellFromXY(mask,coordinates(sp))] <- NA
#
#               # checking of nb candidates
#               nb.cells <- .nb.available.pa.cells(mask)
#               if(nb.cells <= nb.points){
#                 nb.repet <- 1
#                 nb.points <- nb.cells
#                 cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
#               }
#
#               # select cells into raster
#               pa.tab.tmp <- matrix(NA, ncol=nb.repet, nrow=nb.points)
#               for( j in 1:ncol(pa.tab.tmp)){
#                 pa.tab.tmp[,j] <- sampleRandom(x=mask, size=nb.points, cells=T)[,"cell"]
#               }
#
#               # puting cells in good format
#               selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
#               pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
#               colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
#               for( j in 1:ncol(pa.tab)){
#                 pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
#               }
#
#               # puting presences, true absences and pseudo absences together
#               xy <- rbind(coordinates(sp), xyFromCell(mask, selected.cells))
#               sp <- as.numeric(unlist(c(as.vector(sp@data), rep(NA,length(selected.cells))), use.names=FALSE))
#               env <- extract(env, xy)
#
#               pa.tab <- rbind(matrix(TRUE,nrow=(nrow(xy)-length(selected.cells)), ncol=ncol(pa.tab)),
#                               pa.tab)
#
#               return(list(xy = xy,
#                           sp = sp,
#                           env = as.data.frame(env),
#                           pa.tab = as.data.frame(pa.tab)))
#
#             }
#           })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# .arranging.pa.table(pa.data, pa.tab, sp.data=NULL, xy=NULL){
#
#   # transforming sp.data into vector if it's not
#   if(!is.null(sp.data)){ # that means that PA were chosed into explanatories data
#     if(inherits(sp.data, 'SpatialPoints')){
#       xy <- coordinates(sp.data)
#       sp.data <- sp.data@data
#     }
#     if(inherits(sp.data, 'Raster')){
#       xy <- rbind(xyFromCell(sp.data, Which(sp.data >= 1), cells=TRUE), xyFromCell(sp.data, Which(sp.data == 0)))
#       sp.data.tmp <- rep(0,nrow(xy))
#       sp.data.tmp[1:length(Which(sp.data >= 1, cells=TRUE))] <- 1
#       sp.data <- sp.data.tmp
#       rm('sp.data.tmp')
#     }
#   }
#
#   # getting PA selected
#
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# # additional hidden functions
# .allAvailableAbs <- function(data.biomod.species){
#   out <- data.biomod.species
#   if( sum(is.na(out)>0) )
#     out[is.na(out)] <- 0
#   return(out)
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# .is.some.na.in.data <- function(sp){
#   if(is.vector(sp)){
#     if(sum(is.na(sp)) == 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
#
#   if(inherits(sp, 'SpatialPoints')){
#     if(sum(is.na(sp[,1])) == 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
#
#   if(inherits(sp, 'Raster')){
#     if(sp@data@min >= 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
# }
.Rescaler5 <-
function(dataToRescale, ref=NULL, name, original=FALSE, weights=NULL)
{
  compress.arg = ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
#     #preparing data
#     #homogenize the format accross original predictions and new projections 
#     if(!class(dataToRescale)[1]=='RasterLayer'){
#         DataF <- as.data.frame(dataToRescale)         
#         colnames(DataF) <- "DataF"                  
#     } else{
#         names(dataToRescale) <-"DataF"
#         DataF <- stack(dataToRescale) 
#     }
#     
    #Creating or loading the scaling model
  if(original){
      if(! file.exists(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/scaling_models/", sep=""))){
        dir.create(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/scaling_models/", sep=""), showWarnings=F)
      }
#       Rescaling_GLM = glm(ref~DataF, data=DataF, family="binomial", mustart = rep(0.5,length(ref)))
      ## customised wgts
      if(is.null(weights)){
        Prevalence=0.5
        nbPres <- sum(ref, na.rm=TRUE)
        nbAbs <- length(ref) - nbPres
        weights <- rep(1,length(ref))
        
        if(nbAbs > nbPres){ # code absences as 1
          weights[which(ref>0)] <- (Prevalence * nbAbs) / (nbPres * (1-Prevalence))
        } else{ # code presences as 1
          weights[which(ref==0 | is.na(ref))] <- (nbPres * (1-Prevalence)) / (Prevalence * nbAbs)
        }
        
         weights = round(weights[]) # test to remove glm & gam warnings
      }

      

      Rescaling_GLM = glm(ref~pred, data=data.frame(ref=as.numeric(ref), pred = as.numeric(dataToRescale)) , family=binomial(link=probit), x=TRUE, weights=weights)
#       Rescaling_GLM = glm(ref~DataF, data=DataF, family=binomial, weights=wgts)
#       Rescaling_GLM = glm(ref~DataF, data=DataF, family=binomial, weights=wgts)
      eval(parse(text=paste("save(Rescaling_GLM, file='", getwd(),"/",
                            unlist(strsplit(name,'_'))[1], "/models/scaling_models/",
                            name, "_scaled' ,  compress=", compress.arg, ")", sep=""))) 
    } else{
      eval(parse(text=paste("load('", getwd(),"/",unlist(strsplit(name,'_'))[1],
                            "/models/scaling_models/",name,"_scaled')", sep="")))
    }
    #make the scaling prediction
    if(! inherits(dataToRescale, "Raster")){
      RescaledData <- predict(Rescaling_GLM, data.frame(pred=as.numeric(dataToRescale)), type="response")
    } else{
      RescaledData <- predict(dataToRescale, model=Rescaling_GLM, type='response')
    }
    
  
#     cat("\n\t\t original range = ", min(DataF) ," - ", max(DataF), "\t scal ranged = ", min(RescaledData), " - ", max(RescaledData) )

    
# 	  if(class(dataToRescale)[1]=='RasterLayer')  RescaledData <- predict(model=Rescaling_GLM, DataF, type="response")    #rasters
	   
    return(RescaledData)
}


.scaling_model <-
  function(dataToRescale, ref=NULL, ...)
  {
    args <- list(...)
    prevalence <- args$prevalence
    weights <- args$weights
    
#     # if no prevalence define a 0.5 is set.
#     if(is.null(prevalence)) prevalence <- 0.5
    
    # if no weights given, some are created to rise the define prevalence
    if(is.null(weights) & ! is.null(prevalence)){
      nbPres <- sum(ref, na.rm=TRUE)
      nbAbs <- length(ref) - nbPres
      weights <- rep(1,length(ref))
      
      if(nbAbs > nbPres){ # code absences as 1
        weights[which(ref>0)] <- (prevalence * nbAbs) / (nbPres * (1-prevalence))
      } else{ # code presences as 1
        weights[which(ref==0 | is.na(ref))] <- (nbPres * (1-prevalence)) / (prevalence * nbAbs)
      }
      
      weights = round(weights[]) # to remove glm & gam warnings
    } else if(is.null(weights)){
      # only 1 weights vector
      weights <- rep(1,length(ref))
    }
     
    # define a glm to scal prediction from 0 to1 
    scaling_model <- glm(ref~pred, data=data.frame(ref=as.numeric(ref), pred = as.numeric(dataToRescale)) , family=binomial(link=probit), x=TRUE, weights=weights)
    
    return(scaling_model)
  }
##' @name response.plot
##' @title Analysis of the response curves of a model within Biomod
##' @description Depreciated function, please use \code{
##' response.plot2} instead
##' 
##' @param model the model for which you want the response curves to be
##' plotted. Compatible with GAM, GBM, GLM, ANN, CTA, RF, FDA and MARS.
##' @param Data the variables for which you want the response curves to be
##' plotted. A data frame is wanted with one column per variable. They
##' have to have the same names as the ones used to calibrate the model.
##' @param show.variables give in the column numbers of 'Data' for
##' selecting the variables that are wanted for plotting
##' @param save.file can be set to "pdf", "jpeg" or "tiff" to save the
##' plot. Pdf options can be changed by setting the default values of 
##' pdf.options().
##' @param name the name of the file produced if save.file is different to
##' "no" (extensions are already included)
##' @param ImageSize the size of the image in pixels if save.file is
##' different to "no". Affects "jpeg" and "tiff" outputs only. Default if
##' 480 pixels which is the R default. 
##' @param plot if TRUE (the default) then a plot is produced. If not, an
##' array containing predictions is returned (see details)
##' 
##' @details
##' Depreciated function, please use \code{response.plot2}
##' instead.
##' 
##' @author Wilfried Thuiller
##' 
##' @references Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J.
##' R. 2005 The evaluation strip: A new and robust method for plotting
##' predicted responses from species distribution models. Ecological
##' Modelling 186, 280-289.
##' 
##' @keywords plot
##' @keywords models
##' @keywords regression
##' @keywords nonlinear
##' @keywords multivariate
##' @keywords nonparametric
##' @keywords tree
##' 
response.plot <-
  function(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE){

    cat("\n! Deprecated function, please use response.plot2 instead!")
    return(TRUE)
  }


##' @name response.plot2
##' @title Function for for plotting predicted responses from species
##' distribution models in 2 or 3 dimensions
##' 
##' @description Adaptation of the Evaluation Strip method proposed by
##' Elith et al.(2005). This function enables to plot the response curves
##' of a model independently of the algorithm used for building the model.
##' It therefore permits a direct comparison of predicted responses from
##' the different statistical approaches on the same data. 
##'   
##' @param models a character vector specifying the models for which the
##' response curves have to be plotted. Compatible with GAM, GBM, GLM, 
##' ANN, CTA, RF, FDA, MARS and MAXENT.
##' @param Data a dataframe containing the variables for which the
##' response curves have to be plotted. They must have the same names as
##' the ones used to calibrate the model. RasterStack are also supported.
##' @param show.variables the names or the column numbers of 'Data' for
##' selecting the variables to be plotted. By default all columns are
##' taken
##' @param do.bivariate 'logical', if FALSE (default), the predicted
##' response curves are plotted for every singe variable independently 
##' (2 dimension). If TRUE, the predicted response curves are represented
##' in 3 dimentions for all pairs of variables
##' @param fixed.var.metric either 'mean' (default), 'median', 'min' or
##' 'max' specifying the statistic used to fix as constant the remaining
##' variables when the predicted response is estimated for one of the
##' variables
##' @param save.file can be set to "pdf", "jpeg" or "tiff" to save the
##' plot. Pdf options can be changed by setting the default values of 
##' pdf.options().
##' @param name the name of the file produced if save.file is different to
##' "no" (extensions are already included)
##' @param ImageSize the size of the image in pixels if save.file is
##' different to "no". Affects "jpeg" and "tiff" outputs only. Default if
##' 480 pixels which is the R default.
##' @param plot if TRUE (the default) then a plot is produced
##' @param \ldots further arguments :
##'   - \code{data_species} : vector containing data species occurrences.
##'   Have to match with \code{Data}. If given, the statistic used to fix
##'   variables value will be calculated only over presences points.
##'   (Considered only if \code{Data} is under table format)
##'   - \code{col} : vector of colors to be used (only for univariate 
##'   case)
##'   - \code{lty} : vector of lines types to be used
##'   - \code{main} : character, the title of the graph (default one based
##'   on first model class is automatically produced if not referred)
##'   - \code{display_title} : logical, display or not graph title
##'   - \code{legend} : logical, add legend to plot (only for univariate
##'   case)
##' 
##' @details
##' For building the predicted response curves, n-1 variables are set
##' constant to a fixed value (mean, median, min or max i.e 
##' \code{fixed.var.metric} arg) and only the remaining one (resp. 2 for
##' 3D response plot) is varying across its whole range (given by
##' \code{Data}). n the case of categorical variable, the most represented
##' class is taken. The variations observed and the curve thus obtained
##' shows the sensibility of the model to that specific variable. This
##' method does not account for interactions between variables.
##' In the evaluation strip initially proposed by Elith et al. 2005 the
##' remaining variables are set to the mean. 
##' 
##' @return a 4 dimentions array is returned. It contains the necessary
##' outputs to produce the plots. This is useful to make your own custom
##' response plot graphics.
##'   
##' Array returned structure : 
##' 
##' - First dimension: the dimension of the predictions
##' - Second dimension: 2 or 3 columns: The first one (resp. the first 
##' two) is (are) the explanatory variable(s) to plot, the last one, the
##' probability of occurrence
##' - Third dimension: The set of environmental variables for which the
##' response.plot was asked to run.
##' - Fourth dimension:the selected models
##' 
##' @author Wilfried Thuiller, Damien Georges
##' 
##' @references 
##' Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J. R. 2005 The
##' evaluation strip: A new and robust method for plotting predicted
##' responses from species distribution models. Ecological Modelling 186,
##' 280-289.
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}
##' 
##' @keywords plot
##' @keywords models
##' @keywords regression
##' @keywords nonlinear
##' @keywords multivariate
##' @keywords nonparametric
##' @keywords tree
##' 
##' @examples
##' \dontrun{
##' ##' species occurrences
##' DataSpecies <- 
##'   read.csv(
##'     system.file("external/species/mammals_table.csv", package="biomod2"), 
##'     row.names = 1
##'   )
##' head(DataSpecies)
##' ##' the name of studied species
##' myRespName <- 'VulpesVulpes'
##'     
##' ##' the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[, myRespName])
##'     
##' ##' the XY coordinates of species data
##' myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##'
##' myExpl <- 
##'   raster::stack(
##'     system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio12.grd", package = "biomod2")
##'   )
##'
##' ##' 1. Formatting Data
##' myBiomodData <- 
##'   BIOMOD_FormatingData(
##'     resp.var = myResp,
##'     expl.var = myExpl,
##'     resp.xy = myRespXY,
##'     resp.name = myRespName
##'   )
##' 
##' ##' 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' ##' 3. Doing Modelisation
##' myBiomodModelOut <- 
##'   BIOMOD_Modeling(
##'     myBiomodData,
##'     models = c('GLM','RF'),
##'     models.options = myBiomodOption,
##'     NbRunEval = 2,
##'     DataSplit = 80,
##'     VarImport = 0,
##'     models.eval.meth = c('TSS','ROC'),
##'     do.full.models = FALSE,
##'     modeling.id = "test"
##'   )
##' ##' 4. Plot response curves
##' ##' 4.1 Load the models for which we want to extract the predicted
##' ##' response curves
##' myGLMs <- BIOMOD_LoadModels(myBiomodModelOut, models = 'GLM')
##'     
##' ##' 4.2 plot 2D response plots
##' myRespPlot2D <- 
##'   response.plot2(
##'     models = myGLMs,
##'     Data = get_formal_data(myBiomodModelOut, 'expl.var'),
##'     show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
##'     do.bivariate = FALSE,
##'     fixed.var.metric = 'median',
##'     col = c("blue", "red"),
##'     legend = TRUE,
##'     data_species = get_formal_data(myBiomodModelOut, 'resp.var')
##'   )
##'     
##' ##' 4.2 plot 3D response plots
##' ###' here only for a lone model (i.e "VulpesVulpes_PA1_AllData_GLM")
##' myRespPlot3D <- 
##'   response.plot2(
##'   models = myGLMs[1],
##'   Data = get_formal_data(myBiomodModelOut, 'expl.var'), 
##'   show.variables = get_formal_data(myBiomodModelOut, 'expl.var.names'),
##'   do.bivariate = TRUE,
##'   fixed.var.metric = 'median',
##'   data_species = get_formal_data(myBiomodModelOut, 'resp.var'),
##'   display_title = FALSE
##' )
##'     
##' ##' all the values used to produce this plot are stored into the
##' ##' returned object you can redo plots by yourself and customised 
##' ##' them
##' dim(myRespPlot2D)
##' dimnames(myRespPlot2D)
##'     
##' dim(myRespPlot3D)
##' dimnames(myRespPlot3D)
##' }
##' 
response.plot2 <- function(
  models,
  Data,
  show.variables = seq(1:ncol(Data)),
  do.bivariate = FALSE,
  fixed.var.metric = 'mean',
  save.file = "no",
  name = "response_curve",
  ImageSize = 480,
  plot = TRUE,
  ...
){

  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  add.args <- list(...)
  on_0_1000 <- add.args$on_0_1000
  col <- add.args$col
  lty <- add.args$lty
  legend <- add.args$legend
  main <- add.args$main
  display_title <- add.args$display_title
  # par <- add.args$par
  restrict_to_pres_range <- add.args$restrict_to_pres_range
  data_species <- add.args$data_species
  use.formal.names <- add.args$use.formal.names

  if(is.null(on_0_1000)) on_0_1000 <- FALSE
  if(is.null(col)) if(!do.bivariate) col <- "black" else col <- c("red","orange","green")
  if(is.null(lty)) lty <- 1
  if(is.null(legend)) legend <- FALSE
  if(is.null(display_title)) display_title <- TRUE
  if(is.null(restrict_to_pres_range)) restrict_to_pres_range <- FALSE
  if(is.null(use.formal.names)) use.formal.names <- FALSE

  formal_names <- models

  args <- .response.plot2.check.arg(models, Data, show.variables, save.file, name, ImageSize, plot, fixed.var.metric, do.bivariate, add.args)

  models <- args$models
  Data <- args$Data
  show.variables <- args$show.variables
  save.file <- args$save.file
  name <- args$name
  ImageSize <- args$ImageSize
  plot <- args$plot
  fixed.var.metric <- args$fixed.var.metric
  do.bivariate <- args$do.bivariate
  nb.pts <- args$nb.pts

  if(is.null(main)) main <- try(paste("Response curves for ", get(models[1])@resp_name, "'s ",  get(models[1])@model_class,sep=""))
  if(is.null(data_species)) data_species <- rep(1,nrow(Data)) else data_species[data_species!=1 | is.na(data_species)] <- 0


  # 2. build function outputs
  factor_id <- which(sapply(Data,is.factor))

  list.out <- list()

  # Create a ranged data table
  ref_table <- Data[1,,drop=F]
  rownames(ref_table) <- NULL

  for(i in 1:ncol(Data)){
    if(is.numeric(Data[,i])){
      ref_table[,i] <- switch(fixed.var.metric,
                              mean = mean(Data[data_species==1,i]),
                              median = median(Data[data_species==1,i]),
                              min = min(Data[data_species==1,i]),
                              max = max(Data[data_species==1,i]))
    } else{
      # return everytimes the majoritary class
      sum_level <- summary(Data[data_species==1,i], na.rm = TRUE)
      ref_table[,i] <- names(sum_level)[which.max(sum_level)]
    }
  }



  if(plot){
    # X. Open a graphic file for plotting restults
    if(save.file=="pdf") pdf(paste(name, "pdf", sep="."))
    if(save.file=="jpeg") jpeg(paste(name, "jpeg", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="tiff") tiff(paste(name, "tiff", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="postscript") postscript(paste(name, "eps", sep="."))

    # XX. parametrize our plot window

    if(!do.bivariate){
      nb.graphs <- length(show.variables)
    } else{
      nb.graphs <- length(models) *  ( (length(show.variables)-1) * length(show.variables) / 2 )
    }

    if(legend) nb.graphs <- nb.graphs + 1

    if(display_title){
      W.width <- ceiling(sqrt(nb.graphs))
      W.height <- ceiling(nb.graphs/W.width)

      mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE)
      layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))

      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
      polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
      text(x=0.5, y=0.8, pos=1, cex=1.6, labels= main ,col="#4c57eb")
      par(mar = c(2,2,3.5,1))
    } else {
      par(mfrow=.CleverCut(nb.graphs))

    }

  }


  if(!do.bivariate){
    for(vari in show.variables){
      if(plot) {
        if(on_0_1000) ylim <- c(0,1000) else ylim <- c(0,1)

        if(is.factor(Data[,vari])) xlim <- c(1,length(levels(Data[,vari]))) else xlim=c(min(Data[,vari], na.rm=T), max(Data[,vari], na.rm=T))

        plot(0,0,col="white",xlim=xlim, ylim=ylim, main=vari, ann=TRUE, bty="o",xaxs="r", xaxt="s", xaxt=ifelse(is.factor(Data[,vari]),'n','t'),
             xlab="",ylab="")
        rug(Data[ ,vari])

        if(is.factor(Data[,vari])) axis(1, at=seq(xlim[1],xlim[2],1), labels=levels(Data[,vari]))

        # define color vector
        col <- rep(col,length.out=length(models))
        lty <- rep(lty,length.out=length(models))
      }


      # creating Tmp data
      if(is.factor(Data[,vari])) pts.tmp <- as.factor(levels(Data[,vari])) else pts.tmp <- seq(min(Data[,vari]), max(Data[,vari]), length.out=nb.pts)

      Data.r.tmp <- eval(parse(text=paste("cbind(",vari,"=pts.tmp,ref_table[,-which(colnames(ref_table)==vari),drop=F])",sep="")))
      Data.r.tmp <- Data.r.tmp[,colnames(ref_table),drop=F]
      if(length(factor_id)){
        for(f in factor_id){
          Data.r.tmp[,f] <- factor(as.character(Data.r.tmp[,f]), levels=levels(Data[,f]))
        }
      }


      for(model in models){


        # 0. get model
        mod <- get(model)
        mod.name <- ifelse(use.formal.names, formal_names[which(is.element(models, model))], model)

        # cat("\n*** model = ", model, ", mod.name =  ", mod.name)


        # 2. make projections
        proj.tmp <- predict(mod, Data.r.tmp, on_0_1000=on_0_1000, do_check=FALSE)

        # 4. Ploting results
        if(plot ) {
          if(is.factor(Data[,vari])){
            points(pts.tmp[1:length(levels(Data[,vari]))], proj.tmp[1:length(levels(Data[,vari]))], col=col[which(models==model)], lty = lty[which(models==model)])
          } else{
            lines(pts.tmp[1:nb.pts], proj.tmp[1:nb.pts], col=col[which(models==model)], lty = lty[which(models==model)])
          }
        }

        # 5. Storing results
        if(length(list.out[[vari]]) == 0){ #init
          eval(parse(text=paste("list.out[['",vari,"']] <- data.frame(",vari,"=pts.tmp, ",mod.name,"=proj.tmp)",sep="")))
        } else{
          eval(parse(text=paste("list.out[['",vari,"']] <- cbind(list.out[['",vari,"']],",mod.name,"=proj.tmp)",sep="")))
        }

      }

    }
    if(legend & plot){
      plot.new()
      legend(x="center",
             legend = formal_names,
             col = col,
             lty = lty,
             bty = 'n')
    }

  } else{ ## bivariate case
    for(vari1 in show.variables[-length(show.variables)]){
      for(vari2 in show.variables[-(1:which(show.variables == vari1))]){


        # creating Tmp data
        #         if(is.factor(Data[,vari])) pts.tmp <- as.factor(levels(Data[,vari])) else pts.tmp <- seq(min(Data[,vari]), max(Data[,vari]), length.out=nb.pts)

        pts.tmp1 <- rep(seq(min(Data[,vari1]), max(Data[,vari1]), length.out=sqrt(nb.pts)),each=sqrt(nb.pts))
        pts.tmp2 <- rep(seq(min(Data[,vari2]), max(Data[,vari2]), length.out=sqrt(nb.pts)),sqrt(nb.pts))

        Data.r.tmp <- eval(parse(text=paste("cbind(",vari1,"=pts.tmp1,",vari2,"=pts.tmp2, ref_table[,-which(colnames(ref_table)%in% c(vari1,vari2)),drop=F])",sep="")))
        Data.r.tmp <- Data.r.tmp[,colnames(ref_table),drop=F]
        if(length(factor_id)){
          for(f in factor_id){
            Data.r.tmp[,f] <- factor(as.character(Data.r.tmp[,f]), levels=levels(Data[,f]))
          }
        }

        for(model in models){

          # 0. get model
          mod <- get(model)
          mod.name <- ifelse(use.formal.names, formal_names[which(is.element(models, model))], model)

          # 2. make projections
          proj.tmp <- predict(mod, Data.r.tmp, on_0_1000=on_0_1000, do_check=FALSE)

          # 4. Storing results
          vari <- paste(vari1,vari2,sep="_")
          if(length(list.out[[vari]]) == 0){ #init
            eval(parse(text=paste("list.out[['",vari,"']] <- data.frame(",vari1,"=pts.tmp1,",vari2,"=pts.tmp2, ",mod.name,"=proj.tmp)",sep="")))
          } else{
            eval(parse(text=paste("list.out[['",vari,"']] <- cbind(list.out[['",vari,"']],",mod.name,"=proj.tmp)",sep="")))
          }

          # 5. Ploting results
          if(plot) {
            # reformating results to perform a persp plot
            pts.tmp1 <- sort(unique(pts.tmp1))
            pts.tmp2 <- sort(unique(pts.tmp2))
            proj.tmp <- matrix(proj.tmp, ncol=length(pts.tmp2), byrow=FALSE)

            # build color scale
            ncz <- length(pts.tmp2)
            nrz <- length(pts.tmp1)
            # Create a function interpolating colors in the range of specified colors
            jet.colors <- colorRampPalette(col)
            # Generate the desired number of colors from this palette
            nbcol <- 50
            color <- jet.colors(nbcol)
            # Compute the z-value at the facet centres
            zfacet <- proj.tmp[-1, -1] + proj.tmp[-1, -ncz] + proj.tmp[-nrz, -1] + proj.tmp[-nrz, -ncz]
            # Recode facet z-values into color indices
            facetcol <- cut(zfacet, nbcol)

            persp(x=pts.tmp1,y=pts.tmp2,z=proj.tmp, xlab = vari1, ylab=vari2, zlab="pred", zlim=c(0,1), theta = 30, phi = 30,
                  expand = 0.5, col = color[facetcol], ltheta = 120, shade = 0.25, ticktype = "simple", main = formal_names[which(models==model)], cex.main = 0.9, cex.axis=0.7)
          }

        }
      }
    }
  }

  # XXX. Close file
  if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()

  # delete temp files if somes has been created
  if(file.exists(file.path(get(models[1])@resp_name,'RespPlotTmp'))){
    unlink(path.expand(file.path(get(models[1])@resp_name,'RespPlotTmp')), recursive=TRUE, force=TRUE)
  }

  # transform list.out into ggplot firendly shape
  if(do.bivariate){
    gg.out <- .as.ggdat.2D(list.out)
  } else {
    gg.out <- .as.ggdat.1D(list.out)
  }

  invisible(gg.out)
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.response.plot2.check.arg <- function(models, Data, show.variables, save.file, name, ImageSize, plot, fixed.var.metric, do.bivariate, add.args){

  # 1. check add.args
  #   if(sum(! (names(add.args) %in% c("nb.pts","xy"))) > 0){
  #     warning(paste(toString(names(add.args)[which(! (names(add.args) %in% c("nb.pts")))]), " are unknown arguments", sep="" ))
  #   }


  ### check of models args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(!is.character(models)){
    stop("models must be a character vector of models names")
  }

  mod_names <- NULL
  for(mod in models){
    if(!exists(mod)){
      stop("you need to load the selected models!")
    }

    if(!inherits(get(mod), 'biomod2_model')){

      # create a biomod2 modeling object
      mod_tmp <- .Construct.default.biomod2.modeling.obj(get(mod))
      assign(mod_tmp@model_name, mod_tmp, envir = parent.frame(n = 1))
      mod_names <- c(mod_names, mod_tmp@model_name)
    } else{
      mod_names <- c(mod_names, mod)
    }
  }

  models <- mod_names


  ### defining the number split in each variables range =-=-=-=-=- #
  if(!is.null(add.args$nb.pts)){
    if(do.bivariate){
      # total number of points is the square of the difined
      add.args$nb.pts <- add.args$nb.pts^2
    }
  } else{
    if(!do.bivariate){
      add.args$nb.pts <- 100
    } else{
      add.args$nb.pts <- 25^2
    }
  }

  ### check of data args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(inherits(Data,"Raster")){
    cat("\n   > Extracting raster infos..")
    DataTmp <- matrix(0,ncol=nlayers(Data), nrow=add.args$nb.pts)
    colnames(DataTmp) <- names(Data)
    maxVal <- maxValue(Data)
    minVal <- minValue(Data)
    for(i in 1:ncol(DataTmp)){
      DataTmp[,i] <- seq(minVal[i],maxVal[i],length.out=add.args$nb.pts)
    }
    Data <- DataTmp
    rm(list=c('maxVal','minVal','DataTmp'))

  }

  ### check show.variables arg -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if( ( length(show.variables) > ncol(Data) ) | (sum(!(show.variables %in% colnames(Data)))) ) stop("columns wanted in show.variables do not match the data \n")

  # remove factorial var in do.bivariate case
  if(do.bivariate){
    fact_var <- sapply(Data[,show.variables, drop=F], is.factor)
    if(sum(fact_var)>0){
      cat("\n\tFactorial variables have been automatically removed!")
      show.variables <- show.variables[!fact_var]
    }
  }

  ### check save.file arg -=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #






  # TO DO
  return(list(models = models,
              Data = Data,
              show.variables = show.variables,
              save.file = save.file,
              name = name,
              ImageSize = ImageSize,
              plot = plot,
              fixed.var.metric = fixed.var.metric,
              do.bivariate = do.bivariate,
              nb.pts = add.args$nb.pts))
}

###

.Construct.default.biomod2.modeling.obj <- function(mod){

  ## ANN ##
  if(inherits(mod, "nnet")){
    return(new("ANN_biomod2_model",
               model = mod,
               model_name = paste(ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_ANN", sep=""),
               model_class = 'ANN',
               resp_name = ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),
               expl_var_names = ifelse( is.character(attr( mod$terms,"term.labels")), attr( mod$terms,"term.labels"), "") ))
  }


  ## CTA ##
  if(inherits(mod, "rpart")){

    return(new("CTA_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_CTA", sep=""),
               model_class = 'CTA',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## FDA ##
  if(inherits(mod, "fda")){
    return(new("FDA_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_FDA", sep=""),
               model_class = 'FDA',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## GAM ##
  if(inherits(mod, "gam")){
    return(new("GAM_biomod2_model",
               model = mod,
               model_subclass = ifelse(mod$method=="glm.fit","GAM_gam","GAM_mgcv"),
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_GAM", sep=""),
               model_class = 'GAM',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## GBM ##
  if(inherits(mod, "gbm")){
    return(new("GBM_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$Terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_GBM", sep=""),
               model_class = 'GBM',
               resp_name = as.character(mod$Terms[[2]]),
               expl_var_names = attr(mod$Terms,"term.labels")))
  }

  ## GLM ##
  if(inherits(mod, c("glm", "lm"))){
    return(new("GLM_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_GLM", sep=""),
               model_class = 'GLM',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## MARS ##
  if(inherits(mod, "mars")){
    return(new("MARS_biomod2_model",
               model = mod,
               model_name =paste("species_AllData_",as.character(format(Sys.time(), "%OS6")),"_MARS",sep=""),
               model_class = 'MARS',
               resp_name = "species",
               expl_var_names = as.character(colnames(mod$factor))))
  }

  ## RF ##
  if(inherits(mod, "randomForest")){
    return(new("RF_biomod2_model",
               model = mod,
               model_name =paste(ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_RF", sep=""),
               model_class = 'RF',
               resp_name = ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),
               expl_var_names = ifelse( is.character(attr( mod$terms,"term.labels")), attr( mod$terms,"term.labels"), "") ))
  }

  stop("Unknown model class")

}

.as.ggdat.1D <-
  function (rp.dat)
  {
    requireNamespace('dplyr')
    out_ <- 
      bind_rows(
        lapply(
          rp.dat, 
          function(dat_) {
            dat_$id <- rownames(dat_)
            id.col.id <- which(colnames(dat_) == "id")
            expl.dat_ <- dat_ %>% 
              dplyr::select(1, id.col.id) %>%
              tidyr::gather("expl.name", "expl.val", 1)
            pred.dat_ <- dat_ %>% 
              dplyr::select(-1, id.col.id) %>%
              tidyr::gather("pred.name", "pred.val", (1:(ncol(dat_)-2)))
            out.dat_ <- 
              dplyr::full_join(expl.dat_, pred.dat_, by = 'id') %>%
              dplyr::mutate_at(c('expl.name', 'pred.name'), as.character) %>%
              dplyr::mutate_at('expl.val', as.numeric)
            return(out.dat_)
          }
        )
      )
    
    out_ <- 
      out_ %>%
      dplyr::mutate_at('expl.name', factor)
    
    return(out_)
  }


.as.ggdat.2D <- 
  function(rp.dat){
    out_ <- 
      bind_rows(
        lapply(
          rp.dat, 
          function(dat_) {
            dat_$id <- rownames(dat_)
            #   dat_$expl.name <- as.character(dat_$expl.name)
            #   dat_$pred.name <- as.character(dat_$pred.name)
            id.col.id <- which(colnames(dat_) == "id")
            expl1.dat_ <- dat_ %>% 
              dplyr::select(1, id.col.id) %>% 
              tidyr::gather("expl1.name", "expl1.val", 1)
            expl2.dat_ <- dat_ %>% 
              dplyr::select(2, id.col.id) %>% 
              tidyr::gather("expl2.name", "expl2.val", 1)
            pred.dat_ <- dat_ %>% 
              dplyr::select(3, id.col.id) %>% 
              tidyr::gather("pred.name", "pred.val", 1)
            out.dat_  <- 
              dplyr::full_join(
                dplyr::full_join(
                  expl1.dat_, 
                  expl2.dat_,
                  by = 'id'
                ), 
                pred.dat_,
                by = 'id'
              )
            out.dat_ <- out.dat_ %>%
              dplyr::mutate_at(c('expl1.name', 'expl2.name', 'pred.name'), as.character) %>% 
              dplyr::mutate_at(c('expl1.val', 'expl2.val'), as.numeric)
            return(out.dat_)
  }))
  ## ensure that the stips are in the right order
  out_ <- 
    out_ %>%
    dplyr::mutate_at(c('expl1.name', 'expl2.name'), factor)
  return(out_)
}
##' 
##' @name sample.factor.levels
##' @aliases sample.factor.levels
##' 
##' @title Tool to ensure the sampling of all levels of a factorial variable 
##' @description This function will sample randomly an element of each level
##'   of all the factorial variables contains in a Raster* object or a data.frame
##' @author damien g.
##' 
##' @param x         a Raster* object or a data.frame
##' @param mask.out  a Raster/data.frame mask containing area that have already 
##'   been sampled. The factor levels within this mask will not be sampled.
##' @param mask.in   a Raster/list of Raster/data.frame mask (potentially a stack of 
##'   masks) indicating areas were we want to sample our factor level in priority.
##'   Note that if after having explore this masks some levels of the considered
##'   factorial varialble remains unsampled, this levels will be sampled in the 
##'   reference input object (here 'x')
##'      
##' @note 
##'   - The x/mask.out/mask.in should be coherent in term of dimention (same number of 
##'     rows for data.frame and same number of rows, column, identic resolution 
##'     and projection coordinates system for Raster* objects)
##'   - If mask.in contains several masks (RasterStack or multi-column data.frame)
##'     then the order of the mask matter. The mask will be considered successively.
##'     The first will be use prioritarly to sample our variable factor levels and 
##'     so on.
##'   - Raster* masks will be understood as: 
##'       - NA: out of of mask
##'       - not NA: in mask
##'   - data.frame masks will be understood as:
##'       - FALSE: out of mask
##'       - TRUE: in mask
##' 
##' @details In case any factorial variable is found in the input object then 
##'   NULL is returned.
##'   
##' @return a numeric vector the number (cell number for Raster* objects or row 
##'   number for data.frame) where each will refer to a single level of a single
##'   factorial variable.
##'   
##' @examples
##' ## example with raster* object ---------- 
##' library(raster)
##' ## create a factorial raster
##' r1 <- raster()
##' r1[] <- 1; r1[1] <- 2; r1[2:3] <- 3
##' r1 <- as.factor(r1)
##' ## create a continuous raster
##' r2 <- raster()
##' r2[] <- rnorm(ncell(r2))
##' ## pull the raster into a RasterStack
##' stk <- stack(r1, r2)
##' is.factor(stk)
##' 
##' ## define a mask for already sampled points
##' mask.out <- r1
##' mask.out[] <- NA; mask.out[2:3] <- 1
##' 
##' ## define a list of mask where we want to sample in priority
##' mask.in.1 <- mask.in.2 <- r1
##' mask.in.1[1:10] <- NA ## only level 1 should be sampled in this mask
##' mask.in.2[1] <- NA ## only levels 1 and 3 should be sampled in this mask
##' mask.in <- list(mask.in.1 = mask.in.1, 
##'                 mask.in.2 = mask.in.2)
##' 
##' ## test different version of the function
##' sample.factor.levels(stk, mask.out = mask.out)
##' sample.factor.levels(stk, mask.in = mask.in)
##' sample.factor.levels(stk, mask.out = mask.out, mask.in = mask.in)
##' 
sample.factor.levels <- function(x, mask.out = NULL, mask.in = NULL){
  ## make some checking of given parameters
  ## TODO(damien)
  if(inherits(x, 'Raster')){
    fact.level.cells <- .sample.factor.levels.raster(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else if(inherits(x, 'data.frame')){
    fact.level.cells <- .sample.factor.levels.data.frame(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else {
    warning(paste0("\nunsupported input data.",
                   "\nx should be a Raster* object or a data.frame.",
                   "\n NULL returned"))
    return(NULL)
  }
}

.sample.factor.levels.raster <- function(x, mask.out = NULL, mask.in = NULL){
  ## identificate the factorial variables
  fact.var <- which(is.factor(x))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f){
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      ## get the levels of the factor on the full dataset
      fact.level.original <- unlist(raster::levels(subset(x, f)))
      fact.level <- fact.level.original
      cat("\n\t> fact.level for",  names(x)[f], ":\t", paste(fact.level, names(fact.level), sep = ":", collapse = "\t"))
      if(!is.null(mask.out)){ ## mask containing points that have already been sampled
        ## check the levels of the fector that have been already sampled
        fact.levels.sampled <- unlist(levels(as.factor(mask(subset(x, f), mask.out))))
        ## update levels names (lost during mask conversion)
        attr(fact.levels.sampled, "names") <- attr(fact.level.original, "names")[fact.levels.sampled]
        cat("\n\t - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- fact.level[!is.element(fact.level, fact.levels.sampled)]
      }
      if(length(fact.level)){
        ## try first to sample factors in the given masks
        if(!is.null(mask.in)){ ## list of mask we want to sample in (order matter!)
          for(mask.in.id in 1:length(mask.in)){
            ## check that some levels remains to be sampled
            if(length(fact.level)){
              ## update the masked version of the factorial raster
              x.f.masked <- as.factor(mask(subset(x, f), mask.in[[mask.in.id]]))
              x.f.levels <- unlist(levels(x.f.masked))
              ## update levels names (lost during mask conversion)
              attr(x.f.levels, "names") <- attr(fact.level.original, "names")[x.f.levels]              
              ## get the list of levels that coulb be sampled in this mask
              fact.levels.in.m.in <- fact.level[is.element(fact.level, x.f.levels)]
              if(length(fact.levels.in.m.in)){
                cat("\n\t - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  sample(which(x.f.masked[] == fl), 1)
                }))
                ## update the list of factor levels to sample
                fact.level <- fact.level[!is.element(fact.level, fact.levels.in.m.in)]
              }
            } 
          } ## end loop over mask.in
        }
        ## @note if some levels remains unsampled then we will take a random value of
        ## them in the full dataset => !! this should be tricky if mask.in arg is given
        ## because the value will be picked out of mask.in but is necessary to 
        ## ensure that models will run smoothly
        if(length(fact.level)){
          cat("\n\t - levels", fact.level, "will be sampled in the original raster")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl){
            sample(which(subset(x, f)[] == fl), 1)}))
        }
      }
      return(unlist(selected.cells))
    })))
    return(fact.level.cells)
  } else { ## no factorial variable
    return(NULL)
  }
}

.sample.factor.levels.data.frame.old <- function(x){
  ## identificate the factorial variables
  fact.var <- which(sapply(x, is.factor))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(sapply(fact.var, function(f){
      fact.level <- levels(x[, f])
      cat("\n fact.level for",  f, ":", fact.level)
      sapply(fact.level, function(fl){
        sample(which(x[,f] == fl), 1)
      })
    }))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}

.sample.factor.levels.data.frame <- function(x, mask.out = NULL, mask.in = NULL){
  ## identificate the factorial variables
  fact.var <- which(sapply(x, is.factor))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f){
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      ## get the levels of the factor on the full dataset
      fact.level.original <- levels(x[, f])
      fact.level <- fact.level.original
      cat("\n> fact.level for",  colnames(x)[f], ":\t", paste(1:length(fact.level), fact.level, sep = ":", collapse = "\t"))
      if(!is.null(mask.out)){ ## mask containing points that have already been sampled
        ## check the levels of the fector that have been already sampled
        fact.levels.sampled <- unique(na.omit(as.character(x[mask.out[, 1], f])))
        ## remove already sampled points from candidates
        x[mask.out[, 1], ] <- NA
#         ## update levels names (lost during mask conversion)
#         attr(fact.levels.sampled, "names") <- attr(fact.level.original, "names")[fact.levels.sampled]
        cat("\n - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- setdiff(fact.level, fact.levels.sampled)
      }
      if(length(fact.level)){
        ## try first to sample factors in the given masks
        if(!is.null(mask.in)){ ## list of mask we want to sample in (order matter!)
          for(mask.in.id in 1:ncol(mask.in)){
            ## check that some levels remains to be sampled
            if(length(fact.level)){
              ## update the masked version of the factorial raster
              x.f.masked <- as.character(x[, f])
              x.f.masked[!mask.in[, mask.in.id]] <- NA
              x.f.levels <- unique(na.omit(x.f.masked))
              ## get the list of levels that coulb be sampled in this mask
              fact.levels.in.m.in <- intersect(fact.level, x.f.levels)
              if(length(fact.levels.in.m.in)){
                cat("\n - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  candidate.cells <- na.omit(which(x.f.masked[] == fl))
                  selected.cell <- NULL
                  if(length(candidate.cells) == 1){ ## single candiate cell
                    selected.cell <- candidate.cells
                  } else { ## multi candidate cells
                    selected.cell <- sample(candidate.cells, 1)
                  }
                  return(selected.cell)
                }))
                ## update the list of factor levels to sample
                fact.level <- setdiff(fact.level, fact.levels.in.m.in)
              }
            } 
          } ## end loop over mask.in
        }
        ## @note if some levels remains unsampled then we will take a random value of
        ## them in the full dataset => !! this should be tricky if mask.in arg is given
        ## because the value will be picked out of mask.in but is necessary to 
        ## ensure that models will run smoothly
        if(length(fact.level)){
          cat("\n - levels", fact.level, "will be sampled in the original data.frame")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl){
            candidate.cells <- na.omit(which(x[, f] == fl))
            selected.cell <- NULL
            if(length(candidate.cells) <= 1){ ## single candidate cell
              selected.cell <- candidate.cells
            } else { ## multi candidate cells
              selected.cell <- sample(candidate.cells, 1)
            }
            return(selected.cell)
          }))
        }
      }
      return(selected.cells)
    })))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}
##' @name SampleMat2
##' @title Sample binary vector
##' 
##' @description
##' \code{SampleMat2} is an internal \pkg{biomod2} function that can help
##' user to sample a binary vector keeping the same proportion of 0s and 1s
##' than in the initial vector.
##' 
##' @param ref a binary vector
##' @param ratio the proportion of \code{ref} to sample
##' @param as.logi logical, if FALSE (default) id of cell will be return;
##'   if TRUE, logical vector of same length than ref will be return
##'   
##' @details
##' This function can be useful to help users to select a part of initial
##' dataset that will be only kept for all validation procedures. 
##' 
##' @return 
##' A list of 2 elements is returned :
##' 
##' - `calibration` Ids of cells selected for calibration (1st sample)
##' 
##' - `evaluation` Ids of cells selected for evaluation (1st sample
##'   complementary)
##'   
##' @author Damien Georges
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_FormatingData}}
##' @keywords models
##' @keywords formula
##' @keywords options
##' 
##' @export
##' 
##' @examples
##' a <- sample(c(0,1),100, replace=TRUE)
##' SampleMat2(ref=a, ratio=0.7)
##' 
SampleMat2 <- function(
  ref, 
  ratio, 
  as.logi = FALSE
){
  # set a new random seed to ensure that sampling is random (issue when CTA is involved and seed needs to be set to a fix number)
  set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6"))*1000000)
  
  ntot <- length(ref)
  npres<- sum(ref)    
  ncal <- ceiling(ntot*ratio)

  pres <- sample(which(ref==1), ceiling(npres*ratio))
  absc <- sample(which(ref==0), ncal-length(pres))
  
  if(as.logi){
    calib <- rep(FALSE, ntot)
    calib[c(pres,absc)] <- TRUE
    eval <- !calib
  } else{
    calib <- c(pres,absc)
    eval <- (1:ntot)[-c(pres,absc)]
  }
  
  return(list("calibration"=calib, "evaluation"=eval))
}

.scope <-
function(enviroTrain, Smoother, degree)
{
    XXX <- enviroTrain
    deg <- degree
    vnames <- names(XXX[])
    step.list <- as.list(vnames)
    names(step.list) <- vnames
    NbVar <- dim(enviroTrain)[2]
    i <- 1
    while(i <= NbVar) {
        vname <- names(XXX)[i]
        # loops through independent variable names
        junk <- c(paste("1 + ",vname, sep=""))
        # minimum scope
        if(is.numeric(XXX[,i])) {
            junk <- c(junk, paste(Smoother, "(", vname, ",", deg, ")", sep=""))
            junk <- eval(parse(text=paste("~", paste(junk, collapse="+"))))
        }
        else if(is.factor(XXX[,i])) {
            junk <- c(junk, paste(vname, sep=""))
            junk <- eval(parse(text=paste("~", paste(junk, collapse="+"))))
        }
        step.list[[vname]] <- junk
        i <- i + 1
    }
    
    return(step.list)
}

`.scope2` <-
function(enviroTrain, formula, Smoother, degree)
{
  # 0. args checking
  if(is.character(formula)) formula <- as.formula(formula)
  if(!inherits(formula,"formula")) stop("formula must be a formula object")
  
  if(is.matrix(enviroTrain)) enviroTrain <- as.data.frame(enviroTrain)
  
  # 1. detect factoriel variables
  factVar <- as.list(names(enviroTrain))
  factVar <- lapply(factVar, is.factor)
  names(factVar) <- names(enviroTrain)
  
  # 2. create the output squeletom
  step.list <- as.list(attr(terms(formula),"term.labels"))
  
  # 3. filling the output obj
  step.list <- lapply(step.list, function(x){
    junk <- c(paste("~1 + ",x, sep=""))
    if(length(factVar[[x]])){ # x is a simple variable
      if(!factVar[[x]]){ # x is not a factor
        junk <- paste(junk, " + ", Smoother, "(", x, ",", degree, ")", sep="")
      }
    } else{
      junk <- paste(junk, " + ", Smoother, "(", x, ",", degree, ")", sep="")
    }
    return(formula(junk))
  })
  
  names(step.list) <- attr(terms(formula),"term.labels")
  
  return(step.list)
}
.scopeExpSyst <-
function(enviroTrain, mod)
{
    i <- 1
    junk2 <- c()
    while(i <= dim(enviroTrain)[2]) {
        
	      vname <- names(enviroTrain)[i]
	      
        if(mod=="NNET" | mod=="FDA" | mod=="GLMs" | mod=="CTA" | mod=="GBM") junk <- vname
        if(mod == "GLMq") {
            if(is.numeric(enviroTrain[,i]))      junk <- paste(vname, "+I(", vname, "^2)+I(",vname, "^3)", sep="")
            else if(is.factor(enviroTrain[,i]))  junk <- vname
        }
        if(mod == "GLMp") {
            if(is.numeric(enviroTrain[,i]))     junk <- paste(vname, "+I(", vname, "^2)+I(",vname, "^3)+", "poly(", vname, ",2) + poly(", vname, ",3)", sep="")
            else if(is.factor(enviroTrain[,i])) junk <- vname
        }
        junk2 <- c(junk2, junk)
        i <- i + 1
    }

    junk2 <- eval(parse(text=paste("~", paste(junk2, collapse="+"))))
    return(junk2)
}

##' @name sre
##' @title Surface Range Envelope
##' @description 
##' Run a rectilinear surface range envelop (equivalent to BIOCLIM) using
##' the extreme percentiles as recommended by Nix or Busby.
##' The SRE performs a simple analysis of within which range of each
##' variable the data is recorded and renders predictions.
##' 
##' @param Response  a vector (resp. a matrix/data.frame, a 
##' SpatialPointsDataFrame or a RasterLayer ) giving your species'
##' presences and absences data
##' @param Explanatory a matrix (resp. a data.frame, a 
##' SpatialPointsDataFrame or a RasterStack ) containing the environmental
##' variables for the sites given in Response. It must have as many rows
##' as there are elements in Response.
##' @param NewData The data for which you want to render predictions with
##' the sre. It must be a matrix (resp. a data.frame, a 
##' SpatialPointsDataFrame or a RasterStack ) of the same type as the one
##' given in Explanatory and with precisely the same variable names.
##' @param Quant the value defines the most extreme values for each
##' variable not to be taken into account for determining the tolerance
##' boundaries for the considered species.
##' @param return_extremcond boolean, if TRUE a matrix containing extreme
##' conditions supported is returned
##' 
##' @details 
##' The more variables you put in, the more restrictive your model will 
##' be (if non-colinear variables).
##' 
##' This method is very much influenced by the data input, and more
##' specifically by the extremes.
##' 
##' Where a linear model can discriminate the extreme values from the main
##' tendency, the SRE considers it equal as any other data point which
##' leads to notable changes in predictions.
##' 
##' Note that, as a consequence of its functioning, the predictions are
##' directly given in binary, a site being either potentially suitable for
##' all the variables, either out of bounds for at least one variable and
##' therefore considered unsuitable. 
##' 
##' The quants argument determines the threshold at which the data will be
##' taken into account for calibration : the default of 0.05 induces that
##' the 5\% most extreme values will be avoided for each variable on each
##' side of its distribution along the gradient. So it in fact takes 5\%
##' away at each end of the variables distribution, giving a total of 10\%
##' of data not considered.
##' 
##' @return 
##' A vector (resp. a RasterLayer ) of the same length as there are rows
##' in NewData giving the prediction in binary (1=presence, 0=absence)
##' 
##' @author Wilfried Thuiller, Bruno Lafourcade, Damien Georges 
##' @seealso \code{\link[biomod2]{BIOMOD_Modeling}}, 
##' \code{\link[biomod2]{BIOMOD_ModelingOptions}}, 
##' \code{\link[biomod2]{BIOMOD_Projection}}
##' 
##' @keywords models
##' @keywords multivariate
##' 
##' @examples
##' require(raster)
##' ##' species occurrences
##' DataSpecies <- 
##'   read.csv(
##'     system.file("external/species/mammals_table.csv", package = "biomod2"), 
##'     row.names = 1
##'   )
##' head(DataSpecies)
##' 
##' ##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' ##' the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' ##' the XY coordinates of species data
##' myRespXY <- DataSpecies[which(myResp==1),c("X_WGS84","Y_WGS84")]
##' 
##' ##' Environmental variables extracted from BIOCLIM (bio_3, 
##' ##' bio_4, bio_7, bio_11 & bio_12)
##' myExpl <- 
##'   raster::stack(
##'     system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio12.grd", package = "biomod2")
##'   )

##' myResp <- 
##'   raster::reclassify(
##'     subset(myExpl, 1, drop = TRUE), c(-Inf, Inf, 0)
##'   )
##' myResp[cellFromXY(myResp,myRespXY)] <- 1
##' 
##' ##' Compute some SRE for several quantile values
##' sre.100 <- 
##'   sre(
##'     Response = myResp, 
##'     Explanatory = myExpl, 
##'     NewData=myExpl, 
##'     Quant = 0
##'   )
##'   
##' sre.095 <- 
##'   sre(
##'     Response = myResp, 
##'     Explanatory = myExpl, 
##'     NewData=myExpl, 
##'     Quant = 0.025
##'   )
##' 
##' sre.090 <- 
##'   sre(
##'     Response = myResp, 
##'     Explanatory = myExpl, 
##'     NewData=myExpl, 
##'     Quant = 0.05
##'   )
##'   
##' ##' visualise results
##' par(mfrow=c(2,2),mar=c(6, 5, 5, 3))
##' plot(myResp, main = paste(myRespName, "original distrib."))
##' plot(sre.100, main="full data calibration")
##' plot(sre.095, main="95 %")
##' plot(sre.090, main="90 %")
sre <- function(
  Response = NULL, 
  Explanatory = NULL, 
  NewData = NULL, 
  Quant = 0.025, 
  return_extremcond = FALSE
){

  # 1. Checking of input arguments validity
  args <- .check.params.sre(Response, Explanatory, NewData, Quant)

  Response <- args$Response
  Explanatory <- args$Explanatory
  NewData <- args$NewData
  Quant <- args$Quant
  rm("args")


  # 2. Determining suitables conditions and make the projection
  lout <- list()
  if(is.data.frame(Response) | is.matrix(Response)){
    nb.resp <- ncol(Response)
    resp.names <- colnames(Response)
    for(j in 1:nb.resp){
      occ.pts <- which(Response[,j]==1)
      extrem.cond <- t(apply(as.data.frame(Explanatory[occ.pts,]), 2,
                           quantile, probs = c(0 + Quant, 1 - Quant), na.rm = TRUE))

      if(!return_extremcond)
        lout[[j]] <- .sre.projection(NewData, extrem.cond)
    }
  }

  if(inherits(Response, 'Raster')){
    nb.resp <- nlayers(Response)
    resp.names <- names(Response)
    for(j in 1:nb.resp){
      occ.pts <- raster::subset(Response,j, drop=T)
      x.ooc.pts <- Which(occ.pts != 1, cells=TRUE, na.rm=T)
      occ.pts[x.ooc.pts] <- rep(NA, length(x.ooc.pts))
      extrem.cond <- quantile(raster::mask(Explanatory, occ.pts), probs = c(0 + Quant, 1 - Quant), na.rm = TRUE)

      if(!return_extremcond)
        lout[[j]] <- .sre.projection(NewData, extrem.cond)
    }
  }

  if(inherits(Response, 'SpatialPoints')){
    nb.resp <- ncol(Response@data)
    resp.names <- colnames(Response@data)
    for(j in 1:nb.resp){
      occ.pts <- which(Response@data[,j]==1)
      if(is.data.frame(Explanatory) | is.matrix(Explanatory)){
        extrem.cond <- t(apply(as.data.frame(Explanatory[occ.pts,]), 2,
                           quantile, probs = c(0 + Quant, 1 - Quant), na.rm = TRUE))
      } else { if(inherits(Explanatory, 'Raster')){
        maskTmp <- raster::subset(Explanatory,1, drop=T)
#         maskTmp <- reclassify(maskTmp, c(-Inf,Inf,NA))
        maskTmp[] <- NA
        maskTmp[cellFromXY(maskTmp, coordinates(Response)[occ.pts,])] <- 1
        extrem.cond <- quantile(raster::mask(Explanatory, maskTmp), probs = c(0 + Quant, 1 - Quant), na.rm = TRUE)
      } else { if(inherits(Explanatory, 'SpatialPoints')){
        ## May be good to check corespondances of Response and Explanatory variables
        extrem.cond <- t(apply(as.data.frame(Explanatory[occ.pts,]), 2,
                   quantile, probs = c(0 + Quant, 1 - Quant), na.rm = TRUE))
      } else { stop("Unsuported case!") } } }

      if(!return_extremcond)
        lout[[j]] <- .sre.projection(NewData, extrem.cond)
    }
  }

  if(return_extremcond){
    return(as.data.frame(extrem.cond))
  } else{
    # 3. Rearranging the lout object
    if(is.data.frame(NewData)){
      lout <- simplify2array(lout)
      colnames(lout) <- resp.names
    }

    if(inherits(NewData, 'Raster')){
      lout <- stack(lout)
      if(nlayers(lout)==1){
        lout <- raster::subset(lout,1,drop=TRUE)
      }
      names(lout) <- resp.names
    }

    return(lout)
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.check.params.sre <- function(Response = NULL, Explanatory = NULL, NewData = NULL, Quant = 0.025){
  # check quantile validity
  if (Quant >= 0.5 | Quant < 0){
    stop("\n settings in Quant should be a value between 0 and 0.5 ")
  }

  # check compatibility between response and explanatory
  if(is.vector(Response) | is.data.frame(Response) | is.matrix(Response) ){
    Response <- as.data.frame(Response)

    if(!is.vector(Explanatory) & !is.data.frame(Explanatory) & !is.matrix(Explanatory) & !inherits(Explanatory, 'SpatialPoints')){
      stop("If Response variable is a vector, a matrix or a data.frame then Explanatory must also be one")
    } else {
      if(inherits(Explanatory, 'SpatialPoints')){
        Explanatory <- as.data.frame(Explanatory@data)
      }
      Explanatory <- as.data.frame(Explanatory)
      nb.expl.vars <- ncol(Explanatory)
      names.expl.vars <- colnames(Explanatory)
      if(nrow(Response) != nrow(Explanatory)){
        stop("Response and Explanatory variables have not the same number of rows")
      }
    }
  }

  if(inherits(Explanatory, 'SpatialPoints')){
    Explanatory <- as.data.frame(Explanatory@data)
    nb.expl.vars <- ncol(Explanatory)
    names.expl.vars <- colnames(Explanatory)
  }

  if(inherits(Response, 'Raster')){
    if(!inherits(Explanatory, 'Raster')){
      stop("If Response variable is raster object then Explanatory must also be one")
    }
    nb.expl.vars <- nlayers(Explanatory)
    names.expl.vars <- names(Explanatory)
  }

  ## check explanatory variables class
  test_no_factorial_var <- TRUE
  if(is.data.frame(Explanatory)){
    if(any(unlist(lapply(Explanatory, is.factor)))){
      test_no_factorial_var <- FALSE
    }
  } else if (inherits(Explanatory, 'Raster')){
    if(any(is.factor(Explanatory))){
      test_no_factorial_var <- FALSE
    }
  }

  if(!test_no_factorial_var) stop("SRE algorithm does not handle factorial variables")


  # If no NewData given, projection will be done on Explanatory variables
  if(is.null(NewData)){
    NewData <- Explanatory
  } else {
    # check of compatible number of explanatories variable
    if(is.vector(NewData) | is.data.frame(NewData) | is.matrix(NewData)){
      NewData <- as.data.frame(NewData)
      if(sum(!(names.expl.vars %in% colnames(NewData))) > 0 ){
        stop("Explanatory variables names differs in the 2 dataset given")
      }
      NewData <- NewData[,names.expl.vars]
      if(ncol(NewData) != nb.expl.vars){
        stop("Incompatible number of variables in NewData objects")
      }
    } else if(!inherits(NewData, 'Raster')){
      NewData <- stack(NewData)
      if(sum(!(names.expl.vars %in% names(NewData))) > 0 ){
        stop("Explanatory variables names differs in the 2 dataset given")
      }
      NewData <- raster::subset(NewData, names.expl.vars)
      if(nlayers(NewData) != nb.expl.vars ){
        stop("Incompatible number of variables in NewData objects")
      }
    }
  }

  return(list(Response = Response,
              Explanatory = Explanatory,
              NewData = NewData,
              Quant = Quant))

}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.sre.projection <- function(NewData, ExtremCond){

  if(is.data.frame(NewData)|is.matrix(NewData)){
    out <- rep(1,nrow(NewData))
    for(j in 1:ncol(NewData)){
      out <- out * as.numeric(NewData[,j] >= ExtremCond[j,1] &  NewData[,j] <= ExtremCond[j,2])
    }
  }

  if(inherits(NewData, "Raster")){
    out <- reclassify(raster::subset(NewData,1,drop=TRUE), c(-Inf, Inf, 1))
    for(j in 1:nlayers(NewData)){
      out <- out * ( raster::subset(NewData,j,drop=TRUE) >= ExtremCond[j,1] ) * ( raster::subset(NewData,j,drop=TRUE) <= ExtremCond[j,2] )
    }
    out <- raster::subset(out,1,drop=TRUE)
  }

  return(out)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

#
#
# sre <- function (Response = NULL, Explanatory = NULL, NewData = NULL, Quant = 0.025)
# {
#   # check quantile validity
#   if (Quant >= 0.5 | Quant < 0){
#     stop("\n settings in Quant should be a value between 0 and 0.5 ")
#   }
#   quants <- c(0 + Quant, 1 - Quant)
#
#   # puting response in an appropriate form if necessary
#   if(inherits(Response, 'Raster')){
#     stop("Response raster stack not suported yet!")
#   } else {
#     Response <- as.data.frame(Response)
#   }
#
#   if(!inherits(NewData, 'Raster')){
#
#   }
#
#
#   if (class(Explanatory)[1] != "RasterStack") {
#
#     if (is.vector(Explanatory)){
#       Explanatory <- as.data.frame(Explanatory)
#       NewData <- as.data.frame(NewData)
#       names(Explanatory) <- names(NewData) <- "VarTmp"
#     }
#     NbVar <- ncol(Explanatory)
#   }
#       if (class(NewData)[1] != "RasterStack"){
#         Pred <- as.data.frame(matrix(0,
#                                      nr = nrow(NewData),
#                                      nc = ncol(Response),
#                                      dimnames = list(seq(nrow(NewData)), colnames(Response))))
#       }
#
#       for (i in 1:ncol(Response)){
#           ref <- as.data.frame(Explanatory[Response[, i] == 1, ])
#           if(ncol(ref)==1){ names(ref) <- names(Explanatory)}
#           if (class(NewData)[1] == "RasterStack") {
#             # select a lone layer
#             TF <- subset(NewData, 1)
#             # put all cell at 1
#             TF <- TF >= TF@data@min
#           }
#           else TF <- rep(1, nrow(NewData))
#
#           for (j in 1:NbVar) {
#               capQ <- quantile(ref[, j], probs = quants, na.rm = TRUE)
#               if (class(NewData)[1] != "RasterStack") {
#                 TF <- TF * (NewData[, names(ref)[j]] >= capQ[1])
#                 TF <- TF * (NewData[, names(ref)[j]] <= capQ[2])
#               }
#               else {
#                 TFmin <- NewData@layers[[which(NewData@layernames ==
#                   names(ref)[j])]] >= capQ[1]
#                 TFmax <- NewData@layers[[which(NewData@layernames ==
#                   names(ref)[j])]] <= capQ[2]
#                 TF <- TF * TFmin * TFmax
#               }
#           }
#           if (class(TF)[1] != "RasterLayer")
#               Pred[, i] <- TF
#           else Pred <- TF
#       }
#   } else{
#
#   }
#
#   if (class(NewData)[1] != "RasterStack" & ncol(Response) == 1)
#   	Pred <- Pred[[1]]
#
#   return(Pred)
#
# }
#
#
#
.testnull <-
function(object, Prev = 0.5 , dat){
  if( is.finite(object$deviance) & is.finite(object$null.deviance)){
    if(object$deviance != object$null.deviance){
      if(inherits(dat,'Raster')){
        pred <- predict(object = dat, model = object, type = "response")
      } else{
        pred <- predict(object, dat, type = "response")
      }
    }
  }
  
  if(!exists('pred')){
    if(inherits(dat,'Raster')){
      pred <- raster::subset(dat,1,drop=TRUE)
      if(Prev < 0.5) pred <- reclassify(x=pred, rcl=c(-Inf,Inf,0))
      if(Prev >= 0.5) pred <- reclassify(x=pred, rcl=c(-Inf,Inf,1))
    } else{
      if(Prev < 0.5) pred <- rep(0, nrow(dat))
      if(Prev >= 0.5) pred <- rep(1, nrow(dat))      
    }
    
  }
    return(pred)
}

# setGeneric(
#   ".transform.outputs",
#   def =
#     function(
#       modOut,
#       out = 'evaluation',
#       ...
#     ){
#       standardGeneric( ".transform.outputs" )
#     }
# )
#
# setMethod('.transform.outputs', signature(modOut='array'),
#   function(modOut, out = 'evaluation'){
#     # check out attr
#     if(!(out %in% c('evaluation', 'prediction', 'var.import', 'calib.failure', 'models.run', 'prediction.eval' ) )){
#       stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'var.import',
#                                                             'calib.failure', 'models.run', 'prediction.eval' ))))
#     }
#
#     # check dim of input list
#     if(length(dim(modOut)) != 4 ){
#       cat('\n',dim(modOut),'\n')
#       print(dimnames(modOut))
#       warning("Not computed .transform.outputs because of an incompatible input list dimension", immediate=T)
#       return(NULL)
#     }
#
#     if(dim(modOut)[4] == 1 & length(unlist(strsplit(unlist(dimnames(modOut)[4]),'_'))) == 1 ){
#       dataset.names <- 'AllData'
#     } else{
#       if(length(dimnames(modOut)[[4]]) > 0){
#         dataset.names <- unlist(sapply(unlist(dimnames(modOut)[4]), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
#       } else {
#         dataset.names <- paste('PA', 1:dim(modOut)[4])
#       }
#     }
#
#     run.eval.names <- sub('_','',unlist(dimnames(modOut)[3]))
#     mod.names <- unlist(dimnames(modOut)[2])
#
#     if (out=='evaluation'){
#       if( is.null(modOut['evaluation',1,1,1])){ return(NULL) }
#       eval.meth.names <- rownames(as.data.frame(modOut['evaluation',1,1,1]))
#       eval.col.names <- colnames(as.data.frame(modOut['evaluation',1,1,1]))
#
#       eval.out <- array(data = unlist(modOut['evaluation',,,]),
#                         dim = c(length(eval.meth.names),
#                                 length(eval.col.names),
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(eval.meth.names,
#                                          eval.col.names,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(eval.out)
#     }
#
#     if (out=='prediction'){
#       if( is.null(modOut['pred',1,1,1])){ return(NULL) }
#       nb.pts.pred <- length(as.numeric(unlist(modOut['pred',1,1,1])))
#       pred.out <- array(data = unlist(modOut['pred',,,]),
#                         dim = c(nb.pts.pred,
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(NULL,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(pred.out)
#     }
#
#     if (out=='prediction.eval'){
#       if( is.null(modOut['pred.eval',1,1,1])){ return(NULL) }
#       nb.pts.pred.eval <- length(as.numeric(unlist(modOut['pred.eval',1,1,1])))
#       pred.eval.out <- array(data = unlist(modOut['pred.eval',,,]),
#                         dim = c(nb.pts.pred.eval,
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(NULL,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(pred.eval.out)
#     }
#
#     if (out=='var.import'){
#       if( is.null(unlist(modOut['var.import',1,1,1]))){ return(NULL) }
#       nb.var <- length(as.numeric(unlist(modOut['var.import',1,1,1])))
#
#       vi.out <- array(data = unlist(modOut['var.import',,,]),
#                         dim = c(nb.var,
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(paste('Var',1:nb.var,sep=''), # to change
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(vi.out)
#     }
#
#     if (out == 'calib.failure'){
#       cf.out <- unlist(modOut['calib.failure',,,])
#       return(cf.out[!is.null(cf.out)])
#     }
#
#     if (out == 'models.run'){
#       mod.run.out <- unlist(modOut['ModelName',,,])
#       return(mod.run.out[!is.null(mod.run.out)])
#     }
#
#   })
#
# # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# setMethod('.transform.outputs', signature(modOut='list'),
#   function(modOut, out = 'evaluation', dim.names = NULL){
#
#     # check out attr
#     if(!(out %in% c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
#                     'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation' ) )){
#       stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'prediction.eval', 'var.import',
#                                                             'calib.failure', 'models.run', 'EF.prediction',
#                                                             'EF.PCA.median', 'EF.evaluation'))))
#     }
#
#     if(length(modOut) == 1 & length(unlist(strsplit(unlist(names(modOut)),'_'))) == 1 ){
#       dataset.names <- 'AllData'
#     } else{
#       if(is.null(dim.names)){
#         dataset.names <- unlist(sapply(unlist(names(modOut)), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
#       } else{
#         dataset.names <- unlist(dim.names[1])
#       }
#     }
#
#     if(is.null(dim.names)){
#       run.eval.names <- sub('_','',unlist(names(modOut[[1]]))) # may be good here to test that all names are identics
#
#       mod.names <- unlist(names(modOut[[1]][[1]]))
#     } else{
#       run.eval.names <- unlist(dim.names[2])
#       mod.names <- unlist(dim.names[3])
#     }
#
#     if (out=='evaluation'){
#
#       eval.tab <- NULL
#       nb_pa <- length(modOut)
#       nb_run <- length(modOut[[1]])
#       nb_mod <- length(modOut[[1]][[1]])
#
#       for(i in 1:nb_pa){
#         for(j in 1:nb_run){
#           for(k in 1:nb_mod){
#             eval.tab <- modOut[[i]][[j]][[k]][['evaluation']]
#             if(!is.null(eval.tab)){ break }
#           }
#           if(!is.null(eval.tab)){ break }
#         }
#         if(!is.null(eval.tab)){ break }
#       }
#
#       if( is.null(eval.tab)){ return(NULL) }
#
#       eval.meth.names <- rownames(as.data.frame(eval.tab))
#       eval.col.names <- colnames(as.data.frame(eval.tab))
#
#       eval.out <- lapply(names(modOut),function(d1){ # data set
#                     lapply(names(modOut[[d1]]), function(d2){ # run eval
#                       lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#                         if(is.null(modOut[[d1]][[d2]][[d3]][['calib.failure']])){
#                           return(data.frame(modOut[[d1]][[d2]][[d3]][['evaluation']]))
#                         } else { matrix(NA, ncol=length(eval.col.names), nrow=length(eval.meth.names), dimnames=list(eval.meth.names,eval.col.names))}
#                       })
#                     })
#                   })
#
#       eval.out <- array(data = unlist(eval.out),
#                         dim = c(length(eval.meth.names),
#                                 length(eval.col.names),
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(eval.meth.names,
#                                          eval.col.names,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(eval.out)
#     }
#
#     if (out=='prediction'){
#
#       pred.tab <- NULL
#       nb_pa <- length(modOut)
#       nb_run <- length(modOut[[1]])
#       nb_mod <- length(modOut[[1]][[1]])
#
#       for(i in 1:nb_pa){
#         for(j in 1:nb_run){
#           for(k in 1:nb_mod){
#             pred.tab <- modOut[[i]][[j]][[k]][['pred']]
#             if(!is.null(pred.tab)){ break }
#           }
#           if(!is.null(pred.tab)){ break }
#         }
#         if(!is.null(pred.tab)){ break }
#       }
#
#       if( is.null(pred.tab)){ return(NULL) }
#
#
#       nb.pts.pred <- length(as.numeric(pred.tab))
#
#       pred.out <- lapply(names(modOut),function(d1){ # data set
#                     lapply(names(modOut[[d1]]), function(d2){ # run eval
#                       lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#                         if(is.null(modOut[[d1]][[d2]][[d3]][['pred']])){
#                           return(rep(NA,nb.pts.pred))
#                         } else{
#                           return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred']]))
#                         }
#                       })
#                     })
#                   })
#
#       pred.out <- array(data = unlist(pred.out),
#                         dim = c(nb.pts.pred,
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(NULL,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(pred.out)
#     }
#
#     if (out=='prediction.eval'){
#       pred.eval.tab <- NULL
#       nb_pa <- length(modOut)
#       nb_run <- length(modOut[[1]])
#       nb_mod <- length(modOut[[1]][[1]])
#
#       for(i in 1:nb_pa){
#         for(j in 1:nb_run){
#           for(k in 1:nb_mod){
#             pred.eval.tab <- modOut[[i]][[j]][[k]][['pred.eval']]
#             if(!is.null(pred.eval.tab)){ break }
#           }
#           if(!is.null(pred.eval.tab)){ break }
#         }
#         if(!is.null(pred.eval.tab)){ break }
#       }
#
#       if( is.null(pred.eval.tab)){ return(NULL) }
#
#
#       nb.pts.pred.eval <- length(as.numeric(pred.eval.tab))
#
#       pred.eval.out <- lapply(names(modOut),function(d1){ # data set
#                     lapply(names(modOut[[d1]]), function(d2){ # run eval
#                       lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#                         if(is.null(modOut[[d1]][[d2]][[d3]][['pred.eval']])){
#                           return(rep(NA,nb.pts.pred.eval))
#                         } else{
#                           return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred.eval']]))
#                         }
#                       })
#                     })
#                   })
#
#       pred.eval.out <- array(data = unlist(pred.eval.out),
#                         dim = c(nb.pts.pred.eval,
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(NULL,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(pred.eval.out)
#     }
#
#     if (out=='var.import'){
#       vi.tab <- NULL
#       nb_pa <- length(modOut)
#       nb_run <- length(modOut[[1]])
#       nb_mod <- length(modOut[[1]][[1]])
#
#       for(i in 1:nb_pa){
#         for(j in 1:nb_run){
#           for(k in 1:nb_mod){
#             vi.tab <- modOut[[i]][[j]][[k]][['var.import']]
#             if(!is.null(vi.tab)){ break }
#           }
#           if(!is.null(vi.tab)){ break }
#         }
#         if(!is.null(vi.tab)){ break }
#       }
#
#       if( is.null(vi.tab)){ return(NULL) }
#
#       nb.var <- length(as.numeric(unlist(vi.tab)))
#
#       ef.mod <- grep(pattern="EF.",mod.names) # EF models
#       if(length(ef.mod)>0){
#         kept.mod <- mod.names[-ef.mod]
#       } else{
#         kept.mod <- mod.names
#       }
#
#       vi.out <- lapply(names(modOut),function(d1){ # data set
#                   lapply(names(modOut[[d1]]), function(d2){ # run eval
#                     lapply(kept.mod, function(d3){ # models without EF ones
#                       if(is.null(modOut[[d1]][[d2]][[d3]][['var.import']])){
#                         return(rep(NA,nb.var))
#                       } else{
#                         return(as.matrix(modOut[[d1]][[d2]][[d3]][['var.import']]))
#                       }
#                     })
#                   })
#                 })
#
#       vi.out <- array(data = unlist(vi.out),
#                         dim = c(nb.var,
#                                 length(kept.mod),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(names(modOut[[1]][[1]][[1]][['var.import']]), # to change
#                                          kept.mod,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(vi.out)
#     }
#
#     if (out == 'calib.failure'){
#       cf.out <- lapply(names(modOut),function(d1){ # data set
#                   lapply(names(modOut[[d1]]), function(d2){ # run eval
#                     lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#                       return(modOut[[d1]][[d2]][[d3]][['calib.failure']])
#                     })
#                   })
#                 })
#       cf.out <- unlist(cf.out)
#       if(length(cf.out)) cf.out <- na.omit(cf.out)
#       if(length(cf.out)) cf.out <- cf.out[!is.null(cf.out)]
#       if(!length(cf.out)) cf.out <- 'none'
#       return(cf.out)
#     }
#
#     if (out == 'models.run'){
#       mod.run.out <- lapply(names(modOut),function(d1){ # data set
#                   lapply(names(modOut[[d1]]), function(d2){ # run eval
#                     lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#                       return(as.character(modOut[[d1]][[d2]][[d3]][['ModelName']]))
#                     })
#                   })
#                 })
#       mod.run.out <- unlist(mod.run.out)
#       if(length(mod.run.out)) mod.run.out <- na.omit(mod.run.out)
#       if(length(mod.run.out)) mod.run.out <- mod.run.out[!is.null(mod.run.out)]
#       if(!length(mod.run.out)) mod.run.out <- 'none'
#       return(mod.run.out)
#     }
#
#
#     if (out == 'EF.prediction'){
#       if( is.null(modOut[[1]][[1]][[1]][['EM']])){ return(NULL) }
#
#       nb.pts.ef.pred <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['EM']])))
#
#       ef.pred.out <- lapply(1:length(modOut),function(d1){ # data set
#                         lapply(1:length(modOut[[d1]]), function(d2){ # run eval
#                           lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
#                             return(as.numeric(modOut[[d1]][[d2]][[d3]][['EM']]))
#                           })
#                         })
#                       })
#
#       ef.pred.out <- array( data = unlist(ef.pred.out),
#                             dim = c(nb.pts.ef.pred,
#                                     length(modOut[[1]][[1]]),
#                                     length(modOut[[1]]),
#                                     length(modOut)),
#                             dimnames = list(NULL,
#                                              mod.names,
#                                              run.eval.names,
#                                              dataset.names))
#
#      return(ef.pred.out)
#     }
#
#     if (out == 'EF.PCA.median'){
#       if( is.null(modOut[[1]][[1]][[1]][['PCA.median']])){ return(NULL) }
#
#       ef.pca.out <- lapply(1:length(modOut),function(d1){ # data set
#                         lapply(1:length(modOut[[d1]]), function(d2){ # run eval
#                           lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
#                             return(as.character(modOut[[d1]][[d2]][[d3]][['PCA.median']]))
#                           })
#                         })
#                       })
#
#       ef.pca.out <- array( data = unlist(ef.pca.out),
#                             dim = c(1,
#                                     length(modOut[[1]][[1]]),
#                                     length(modOut[[1]]),
#                                     length(modOut)),
#                             dimnames = list(NULL,
#                                              mod.names,
#                                              run.eval.names,
#                                              dataset.names))
#
#      return(ef.pca.out)
#     }
#
#     if (out == 'EF.evaluation'){
#       if( is.null(modOut[[1]][[1]][[1]][['EM.eval']])){ return(NULL) }
#       eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
#       eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
#
#       ef.eval.out <- lapply(1:length(modOut),function(d1){ # data set
#                     lapply(1:length(modOut[[d1]]), function(d2){ # run eval
#                       lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
#                         return(data.frame(modOut[[d1]][[d2]][[d3]][['EM.eval']]))
#                       })
#                     })
#                   })
#
#       ef.eval.out <- array(data = unlist(ef.eval.out),
#                         dim = c(length(eval.meth.names),
#                                 length(eval.col.names),
#                                 length(modOut[[1]][[1]]),
#                                 length(modOut[[1]]),
#                                 length(modOut)),
#                         dimnames = list(eval.meth.names,
#                                          eval.col.names,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
#
#      return(ef.eval.out)
#     }
#
#   })
#
# DF_to_ARRAY <- function(df){
# #   cat("\n*** class(df) = ", class(df))
# #   cat("\n*** colnames(df) = ", colnames(df))
#   if(!is.data.frame(df) & !is.matrix(df)){
#     if(is.list(df)){
#       df.names <- names(df)
#       df <- as.data.frame(df)
#       names(df) <- df.names
#     } else{
#       stop("You have to give a data.frame")
#     }
#   }
#
#   a <- sapply(strsplit(colnames(df), '_'), tail, n=3)
#   b <- lapply(1:3, function(id) return(unique(a[id,])))
#   array.dim.names <- c(list(character(0)),rev(b))
# #   array.dim.names <- c(list(c(NULL)),rev(apply(sapply(strsplit(colnames(df), '_'), tail, n=3),1,unique)))
#
#   array.dim <- c(nrow(df),sapply(array.dim.names[-1],length))
#   array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
#
#   for(x in colnames(df)){
#     dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
#     array.out[,dimTmp[1],dimTmp[2],dimTmp[3]] <- df[,x]
#   }
#   return(array.out)
# }
#
# LIST_to_ARRAY <- function(ll){
#   test <- sapply(ll, is.array)
#   if(!all(test)) stop("list elements should be arrays")
#   test <- sapply(ll,dim)
#   test <- apply(test,1,function(x){length(unique(x))==1})
#   if(!all(test)) stop("list elements differ in dimension")
#
#   formal.dim.names <- dimnames(ll[[1]])
#   new.dim.names <- rev(apply(sapply(strsplit(names(ll), '_'), tail, n=3),1,unique))
#   array.dim.names <- c(formal.dim.names,new.dim.names)
#   array.dim <- sapply(array.dim.names,length)
#
#   array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
#
#   for(x in names(ll)){
#     dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
#     dimTmp <- paste( paste(rep(",",length(formal.dim.names)),collapse="") , paste("'",dimTmp,"'",sep="",collapse=","),collapse="")
#
#     eval(parse(text=paste("array.out[",dimTmp,"] <-  ll[[x]]",sep="")))
#   }
#   return(array.out)
# }


# #' Reshape biomod2 objects
# #'
# #' This is an internal function (developper only)
# #'
# #' @inheritParams .transform.outputs.list
# #'
# #' @return extracted statistics of interest from the model object
# #'   as `array`.
# #' @export
# #'

# .transform.outputs.array <-
#   function(
#     modOut,
#     out = 'evaluation'
#   ){
#     # check out attr
#     if(!(out %in% c('evaluation', 'prediction', 'var.import', 'calib.failure', 'models.run', 'prediction.eval' ) )){
#       stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'var.import',
#                                                             'calib.failure', 'models.run', 'prediction.eval' ))))
#     }
# 
#     # check dim of input list
#     if(length(dim(modOut)) != 4 ){
#       cat('\n',dim(modOut),'\n')
#       print(dimnames(modOut))
#       warning("Not computed .transform.outputs because of an incompatible input list dimension", immediate=T)
#       return(NULL)
#     }
# 
#     if(dim(modOut)[4] == 1 & length(unlist(strsplit(unlist(dimnames(modOut)[4]),'_'))) == 1 ){
#       dataset.names <- 'AllData'
#     } else{
#       if(length(dimnames(modOut)[[4]]) > 0){
#         dataset.names <- unlist(sapply(unlist(dimnames(modOut)[4]), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
#       } else {
#         dataset.names <- paste('PA', 1:dim(modOut)[4])
#       }
#     }
# 
#     run.eval.names <- sub('_','',unlist(dimnames(modOut)[3]))
#     mod.names <- unlist(dimnames(modOut)[2])
# 
#     if (out=='evaluation'){
#       if( is.null(modOut['evaluation',1,1,1])){ return(NULL) }
#       eval.meth.names <- rownames(as.data.frame(modOut['evaluation',1,1,1]))
#       eval.col.names <- colnames(as.data.frame(modOut['evaluation',1,1,1]))
# 
#       eval.out <- array(data = unlist(modOut['evaluation',,,]),
#                         dim = c(length(eval.meth.names),
#                                 length(eval.col.names),
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(eval.meth.names,
#                                         eval.col.names,
#                                         mod.names,
#                                         run.eval.names,
#                                         dataset.names))
# 
#       return(eval.out)
#     }
# 
#     if (out=='prediction'){
#       if( is.null(modOut['pred',1,1,1])){ return(NULL) }
#       nb.pts.pred <- length(as.numeric(unlist(modOut['pred',1,1,1])))
#       pred.out <- array(data = unlist(modOut['pred',,,]),
#                         dim = c(nb.pts.pred,
#                                 length(mod.names),
#                                 length(run.eval.names),
#                                 length(dataset.names)),
#                         dimnames = list(NULL,
#                                         mod.names,
#                                         run.eval.names,
#                                         dataset.names))
# 
#       return(pred.out)
#     }
# 
#     if (out=='prediction.eval'){
#       if( is.null(modOut['pred.eval',1,1,1])){ return(NULL) }
#       nb.pts.pred.eval <- length(as.numeric(unlist(modOut['pred.eval',1,1,1])))
#       pred.eval.out <- array(data = unlist(modOut['pred.eval',,,]),
#                              dim = c(nb.pts.pred.eval,
#                                      length(mod.names),
#                                      length(run.eval.names),
#                                      length(dataset.names)),
#                              dimnames = list(NULL,
#                                              mod.names,
#                                              run.eval.names,
#                                              dataset.names))
# 
#       return(pred.eval.out)
#     }
# 
#     if (out=='var.import'){
#       if( is.null(unlist(modOut['var.import',1,1,1]))){ return(NULL) }
#       nb.var <- length(as.numeric(unlist(modOut['var.import',1,1,1])))
# 
#       vi.out <- array(data = unlist(modOut['var.import',,,]),
#                       dim = c(nb.var,
#                               length(mod.names),
#                               length(run.eval.names),
#                               length(dataset.names)),
#                       dimnames = list(paste('Var',1:nb.var,sep=''), # to change
#                                       mod.names,
#                                       run.eval.names,
#                                       dataset.names))
# 
#       return(vi.out)
#     }
# 
#     if (out == 'calib.failure'){
#       cf.out <- unlist(modOut['calib.failure',,,])
#       return(cf.out[!is.null(cf.out)])
#     }
# 
#     if (out == 'models.run'){
#       mod.run.out <- unlist(modOut['ModelName',,,])
#       return(mod.run.out[!is.null(mod.run.out)])
#     }
# 
#   }

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

#
# #' Reshape biomod2 objects
# #'
# #' This is an internal function (developper only)
# #'
# #' @param modOut the object to transform given as a list
# #' @param out character, the type of input object
# #' @param dim.names character, if not `NULL` the resshaped object will be stored on the hard drive
# #'
# #' @return
# #' @export
# #'
# #' @examples
# 
# .transform.outputs.list =
#   function(
#     modOut,
#     out = 'evaluation',
#     dim.names = NULL
#   ){
# 
#   # check out attr
#   if(!(out %in% c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
#                   'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation' ) )){
#     stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'prediction.eval', 'var.import',
#                                                           'calib.failure', 'models.run', 'EF.prediction',
#                                                           'EF.PCA.median', 'EF.evaluation'))))
#   }
# 
#   if(length(modOut) == 1 & length(unlist(strsplit(unlist(names(modOut)),'_'))) == 1 ){
#     dataset.names <- 'AllData'
#   } else{
#     if(is.null(dim.names)){
#       dataset.names <- unlist(sapply(unlist(names(modOut)), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
#     } else{
#       dataset.names <- unlist(dim.names[1])
#     }
#   }
# 
#   if(is.null(dim.names)){
#     run.eval.names <- sub('_','',unlist(names(modOut[[1]]))) # may be good here to test that all names are identics
# 
#     mod.names <- unlist(names(modOut[[1]][[1]]))
#   } else{
#     run.eval.names <- unlist(dim.names[2])
#     mod.names <- unlist(dim.names[3])
#   }
# 
#   if (out=='evaluation'){
# 
#     eval.tab <- NULL
#     nb_pa <- length(modOut)
#     nb_run <- length(modOut[[1]])
#     nb_mod <- length(modOut[[1]][[1]])
# 
#     for(i in 1:nb_pa){
#       for(j in 1:nb_run){
#         for(k in 1:nb_mod){
#           eval.tab <- modOut[[i]][[j]][[k]][['evaluation']]
#           if(!is.null(eval.tab)){ break }
#         }
#         if(!is.null(eval.tab)){ break }
#       }
#       if(!is.null(eval.tab)){ break }
#     }
# 
#     if( is.null(eval.tab)){ return(NULL) }
# 
#     eval.meth.names <- rownames(as.data.frame(eval.tab))
#     eval.col.names <- colnames(as.data.frame(eval.tab))
# 
#     eval.out <- lapply(names(modOut),function(d1){ # data set
#       lapply(names(modOut[[d1]]), function(d2){ # run eval
#         lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#           if(is.null(modOut[[d1]][[d2]][[d3]][['calib.failure']])){
#             return(data.frame(modOut[[d1]][[d2]][[d3]][['evaluation']]))
#           } else { matrix(NA, ncol=length(eval.col.names), nrow=length(eval.meth.names), dimnames=list(eval.meth.names,eval.col.names))}
#         })
#       })
#     })
# 
#     eval.out <- array(data = unlist(eval.out),
#                       dim = c(length(eval.meth.names),
#                               length(eval.col.names),
#                               length(mod.names),
#                               length(run.eval.names),
#                               length(dataset.names)),
#                       dimnames = list(eval.meth.names,
#                                       eval.col.names,
#                                       mod.names,
#                                       run.eval.names,
#                                       dataset.names))
# 
#     return(eval.out)
#   }
# 
#   if (out=='prediction'){
# 
#     pred.tab <- NULL
#     nb_pa <- length(modOut)
#     nb_run <- length(modOut[[1]])
#     nb_mod <- length(modOut[[1]][[1]])
# 
#     for(i in 1:nb_pa){
#       for(j in 1:nb_run){
#         for(k in 1:nb_mod){
#           pred.tab <- modOut[[i]][[j]][[k]][['pred']]
#           if(!is.null(pred.tab)){ break }
#         }
#         if(!is.null(pred.tab)){ break }
#       }
#       if(!is.null(pred.tab)){ break }
#     }
# 
#     if( is.null(pred.tab)){ return(NULL) }
# 
# 
#     nb.pts.pred <- length(as.numeric(pred.tab))
# 
#     pred.out <- lapply(names(modOut),function(d1){ # data set
#       lapply(names(modOut[[d1]]), function(d2){ # run eval
#         lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#           if(is.null(modOut[[d1]][[d2]][[d3]][['pred']])){
#             return(rep(NA,nb.pts.pred))
#           } else{
#             return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred']]))
#           }
#         })
#       })
#     })
# 
#     pred.out <- array(data = unlist(pred.out),
#                       dim = c(nb.pts.pred,
#                               length(mod.names),
#                               length(run.eval.names),
#                               length(dataset.names)),
#                       dimnames = list(NULL,
#                                       mod.names,
#                                       run.eval.names,
#                                       dataset.names))
# 
#     return(pred.out)
#   }
# 
#   if (out=='prediction.eval'){
#     pred.eval.tab <- NULL
#     nb_pa <- length(modOut)
#     nb_run <- length(modOut[[1]])
#     nb_mod <- length(modOut[[1]][[1]])
# 
#     for(i in 1:nb_pa){
#       for(j in 1:nb_run){
#         for(k in 1:nb_mod){
#           pred.eval.tab <- modOut[[i]][[j]][[k]][['pred.eval']]
#           if(!is.null(pred.eval.tab)){ break }
#         }
#         if(!is.null(pred.eval.tab)){ break }
#       }
#       if(!is.null(pred.eval.tab)){ break }
#     }
# 
#     if( is.null(pred.eval.tab)){ return(NULL) }
# 
# 
#     nb.pts.pred.eval <- length(as.numeric(pred.eval.tab))
# 
#     pred.eval.out <- lapply(names(modOut),function(d1){ # data set
#       lapply(names(modOut[[d1]]), function(d2){ # run eval
#         lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#           if(is.null(modOut[[d1]][[d2]][[d3]][['pred.eval']])){
#             return(rep(NA,nb.pts.pred.eval))
#           } else{
#             return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred.eval']]))
#           }
#         })
#       })
#     })
# 
#     pred.eval.out <- array(data = unlist(pred.eval.out),
#                            dim = c(nb.pts.pred.eval,
#                                    length(mod.names),
#                                    length(run.eval.names),
#                                    length(dataset.names)),
#                            dimnames = list(NULL,
#                                            mod.names,
#                                            run.eval.names,
#                                            dataset.names))
# 
#     return(pred.eval.out)
#   }
# 
#   if (out=='var.import'){
#     vi.tab <- NULL
#     nb_pa <- length(modOut)
#     nb_run <- length(modOut[[1]])
#     nb_mod <- length(modOut[[1]][[1]])
# 
#     for(i in 1:nb_pa){
#       for(j in 1:nb_run){
#         for(k in 1:nb_mod){
#           vi.tab <- modOut[[i]][[j]][[k]][['var.import']]
#           if(!is.null(vi.tab)){ break }
#         }
#         if(!is.null(vi.tab)){ break }
#       }
#       if(!is.null(vi.tab)){ break }
#     }
# 
#     if( is.null(vi.tab)){ return(NULL) }
# 
#     nb.var <- length(as.numeric(unlist(vi.tab)))
# 
#     ef.mod <- grep(pattern="EF.",mod.names) # EF models
#     if(length(ef.mod)>0){
#       kept.mod <- mod.names[-ef.mod]
#     } else{
#       kept.mod <- mod.names
#     }
# 
#     vi.out <- lapply(names(modOut),function(d1){ # data set
#       lapply(names(modOut[[d1]]), function(d2){ # run eval
#         lapply(kept.mod, function(d3){ # models without EF ones
#           if(is.null(modOut[[d1]][[d2]][[d3]][['var.import']])){
#             return(rep(NA,nb.var))
#           } else{
#             return(as.matrix(modOut[[d1]][[d2]][[d3]][['var.import']]))
#           }
#         })
#       })
#     })
# 
#     vi.out <- array(data = unlist(vi.out),
#                     dim = c(nb.var,
#                             length(kept.mod),
#                             length(run.eval.names),
#                             length(dataset.names)),
#                     dimnames = list(names(modOut[[1]][[1]][[1]][['var.import']]), # to change
#                                     kept.mod,
#                                     run.eval.names,
#                                     dataset.names))
# 
#     return(vi.out)
#   }
# 
#   if (out == 'calib.failure'){
#     cf.out <- lapply(names(modOut),function(d1){ # data set
#       lapply(names(modOut[[d1]]), function(d2){ # run eval
#         lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#           return(modOut[[d1]][[d2]][[d3]][['calib.failure']])
#         })
#       })
#     })
#     cf.out <- unlist(cf.out)
#     if(length(cf.out)) cf.out <- na.omit(cf.out)
#     if(length(cf.out)) cf.out <- cf.out[!is.null(cf.out)]
#     if(!length(cf.out)) cf.out <- 'none'
#     return(cf.out)
#   }
# 
#   if (out == 'models.run'){
#     mod.run.out <- lapply(names(modOut),function(d1){ # data set
#       lapply(names(modOut[[d1]]), function(d2){ # run eval
#         lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
#           return(as.character(modOut[[d1]][[d2]][[d3]][['ModelName']]))
#         })
#       })
#     })
#     mod.run.out <- unlist(mod.run.out)
#     if(length(mod.run.out)) mod.run.out <- na.omit(mod.run.out)
#     if(length(mod.run.out)) mod.run.out <- mod.run.out[!is.null(mod.run.out)]
#     if(!length(mod.run.out)) mod.run.out <- 'none'
#     return(mod.run.out)
#   }
# 
# 
#   if (out == 'EF.prediction'){
#     if( is.null(modOut[[1]][[1]][[1]][['EM']])){ return(NULL) }
# 
#     nb.pts.ef.pred <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['EM']])))
# 
#     ef.pred.out <- lapply(1:length(modOut),function(d1){ # data set
#       lapply(1:length(modOut[[d1]]), function(d2){ # run eval
#         lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
#           return(as.numeric(modOut[[d1]][[d2]][[d3]][['EM']]))
#         })
#       })
#     })
# 
#     ef.pred.out <- array( data = unlist(ef.pred.out),
#                           dim = c(nb.pts.ef.pred,
#                                   length(modOut[[1]][[1]]),
#                                   length(modOut[[1]]),
#                                   length(modOut)),
#                           dimnames = list(NULL,
#                                           mod.names,
#                                           run.eval.names,
#                                           dataset.names))
# 
#     return(ef.pred.out)
#   }
# 
#   if (out == 'EF.PCA.median'){
#     if( is.null(modOut[[1]][[1]][[1]][['PCA.median']])){ return(NULL) }
# 
#     ef.pca.out <- lapply(1:length(modOut),function(d1){ # data set
#       lapply(1:length(modOut[[d1]]), function(d2){ # run eval
#         lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
#           return(as.character(modOut[[d1]][[d2]][[d3]][['PCA.median']]))
#         })
#       })
#     })
# 
#     ef.pca.out <- array( data = unlist(ef.pca.out),
#                          dim = c(1,
#                                  length(modOut[[1]][[1]]),
#                                  length(modOut[[1]]),
#                                  length(modOut)),
#                          dimnames = list(NULL,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
# 
#     return(ef.pca.out)
#   }
# 
#   if (out == 'EF.evaluation'){
#     if( is.null(modOut[[1]][[1]][[1]][['EM.eval']])){ return(NULL) }
#     eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
#     eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
# 
#     ef.eval.out <- lapply(1:length(modOut),function(d1){ # data set
#       lapply(1:length(modOut[[d1]]), function(d2){ # run eval
#         lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
#           return(data.frame(modOut[[d1]][[d2]][[d3]][['EM.eval']]))
#         })
#       })
#     })
# 
#     ef.eval.out <- array(data = unlist(ef.eval.out),
#                          dim = c(length(eval.meth.names),
#                                  length(eval.col.names),
#                                  length(modOut[[1]][[1]]),
#                                  length(modOut[[1]]),
#                                  length(modOut)),
#                          dimnames = list(eval.meth.names,
#                                          eval.col.names,
#                                          mod.names,
#                                          run.eval.names,
#                                          dataset.names))
# 
#     return(ef.eval.out)
#   }
# 
# }

# DF_to_ARRAY <- function(df){
#   #   cat("\n*** class(df) = ", class(df))
#   #   cat("\n*** colnames(df) = ", colnames(df))
#   if(!is.data.frame(df) & !is.matrix(df)){
#     if(is.list(df)){
#       df.names <- names(df)
#       df <- as.data.frame(df)
#       names(df) <- df.names
#     } else{
#       stop("You have to give a data.frame")
#     }
#   }
# 
#   a <- sapply(strsplit(colnames(df), '_'), tail, n=3)
#   b <- lapply(1:3, function(id) return(unique(a[id,])))
#   array.dim.names <- c(list(character(0)),rev(b))
#   #   array.dim.names <- c(list(c(NULL)),rev(apply(sapply(strsplit(colnames(df), '_'), tail, n=3),1,unique)))
# 
#   array.dim <- c(nrow(df),sapply(array.dim.names[-1],length))
#   array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
# 
#   for(x in colnames(df)){
#     dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
#     array.out[,dimTmp[1],dimTmp[2],dimTmp[3]] <- df[,x]
#   }
#   return(array.out)
# }
# 
# LIST_to_ARRAY <- function(ll){
#   test <- sapply(ll, is.array)
#   if(!all(test)) stop("list elements should be arrays")
#   test <- sapply(ll,dim)
#   test <- apply(test,1,function(x){length(unique(x))==1})
#   if(!all(test)) stop("list elements differ in dimension")
# 
#   formal.dim.names <- dimnames(ll[[1]])
#   new.dim.names <- rev(apply(sapply(strsplit(names(ll), '_'), tail, n=3),1,unique))
#   array.dim.names <- c(formal.dim.names,new.dim.names)
#   array.dim <- sapply(array.dim.names,length)
# 
#   array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
# 
#   for(x in names(ll)){
#     dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
#     dimTmp <- paste( paste(rep(",",length(formal.dim.names)),collapse="") , paste("'",dimTmp,"'",sep="",collapse=","),collapse="")
# 
#     eval(parse(text=paste("array.out[",dimTmp,"] <-  ll[[x]]",sep="")))
#   }
#   return(array.out)
# }
# 
# 
##' @name update_objects
##' @aliases  update_objects
##' @title update biomod2 objects
##' @description 
##' This function is wrapper to update objects construct with a old
##' version of biomod2 to a current one
##' 
##' @param obj a \code{biomod2} object you want to update (e.g. 
##'   'BIOMOD.formated.data', 'biomod2_model' object)
##' @param recursive logical, if TRUE all objects on which updated object is based will be also updated
##' 
##' @details 
##' This function will add/change/delete all object slots that have
##' evolved between 2 versions of the package.
##' If required, objects stored on hard drive will also be updated.
##' 
##' @return the updated version of the  biomod2 object is return
##' 
##' @author Damien Georges
##' @seealso \code{\link{variables_importance}}, 
##' \code{\link{full_suffling}}
##' 
update_objects <- function(obj, recursive=TRUE){
  if(inherits(obj,'BIOMOD.formated.data') | inherits(obj,'BIOMOD.formated.data.PA')){
    cat("\n=-=-=- BIOMOD.formated.data conversion")
    update.objects_BIOMOD.formated.data(obj, recursive=recursive)
  }

  if(inherits(obj,'BIOMOD.models.out')){
    cat("\n=-=-=- BIOMOD.models.out conversion")
    update.objects_BIOMOD.models.out(obj, recursive=recursive)
  }

}

update.objects_BIOMOD.formated.data <- function(obj,recursive=TRUE){
  slots_to_add <- test_slots(obj)
  if(length(slots_to_add)){
    ref_obj <- new(class(obj))
    slot_to_add_str <- paste(slots_to_add,"=ref_obj@",slots_to_add, sep="", collapse=",")
    eval(parse(text= paste ( "obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     eval(parse(text= paste ( "current_obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     obj <- current_obj
  } else {
    cat("\tnothing to do!\n")
  }
  return(obj)
}

update.objects_biomod2_model <- function(obj,model.out=NULL){
  slots_to_add <- test_slots(obj)
  if(length(slots_to_add)){
    ref_obj <- new(class(obj))

    ## special args
    if(sum(c("expl_var_names","expl_var_type", "expl_var_range") %in% slots_to_add, na.rm=T) > 0){
      if(length(model.out)){
        data <- getModelsInputData(model.out,'expl.var')

        if("expl_var_names" %in% slots_to_add) ref_obj@expl_var_names <- names(data)
        if("expl_var_type" %in% slots_to_add) ref_obj@expl_var_type <- get_var_type(data)
        if("expl_var_range" %in% slots_to_add) ref_obj@expl_var_range <- get_var_range(data)
      }
    }

    slot_to_add_str <- paste(slots_to_add,"=ref_obj@",slots_to_add, sep="", collapse=",")
    eval(parse(text= paste ( "obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     eval(parse(text= paste ( "current_obj <- new('",class(obj),"', obj,",slot_to_add_str,")",sep="")))
#     obj <- current_obj
  } else {
    cat("\n\tnothing to do!\n")
  }
  return(obj)
}

update.objects_BIOMOD.models.out <- function(obj, recursive=T){
  slots_to_add <- test_slots(obj)
  if(length(slots_to_add)){
    ref_obj <- new("BIOMOD.models.out")
    slot_to_add_str <- paste(slots_to_add,"=ref_obj@",slots_to_add, sep="", collapse=",")
    eval(parse(text= paste ( "current_obj <- new('BIOMOD.models.out', obj,",slot_to_add_str,")",sep="")))
  }

  if(recursive){
    cat("\n\tBIOMOD.formated.data checking...")
    data <- getModelsInputData(obj)
    data <- update.objects_BIOMOD.formated.data(data)
    save(data, file=obj@formated.input.data@link)

    cat("\n\tIndividual models checking...\n")
    mod_to_check <- BIOMOD_LoadModels(obj)
    for(mtc in mod_to_check){
      cat("\t",mtc,"\n")
      assign(mtc,update.objects_biomod2_model(get(mtc),model.out=obj))
      save(list=mtc,file=file.path(obj@sp.name,"models",obj@modeling.id, mtc))
    }
  }

  obj_name <- tail(unlist(strsplit(obj@link,.Platform$file.sep,fixed=T)),1)
  assign(obj_name, obj)
  save( list=obj_name,file=obj@link)

}

test_slots <- function(obj){
  out <- NULL
  for(slot in slotNames(obj)){
    test <- try(getElement(obj,slot),silent=T)
    if(inherits(test,"try-error")){
      cat("\t\tredefining",slot,"slot\n")
      out <- c(out, slot)
    }
  }
  return(out)
}
##' @name variables_importance
##' @aliases variables_importance
##' @title Variables importance calculation
##' @description
##' This function will return a variable importance value for 
##' each variable involved within your model.
##' 
##' @param model the model you want to study variables importance
##'   (one of the models supported within biomod2, ensemble models
##'    are also supported
##' @param data the \code{data.set} on which you want to perform
##'   analyses
##' @param method the randomisation method (only 'full_rand'
##'   available so far)
##' @param nb_rand the number of permutation done for each
##'   variable
##' @param ... additional args (not implemented yet)
##' 
##' @details
##' It's more or less base on the same principle than 
##' \code{\link[randomForest]{randomForest}} variables importance
##' algorithm. The principle is to shuffle a single variable of
##' the given data. Make model prediction with this 'shuffled'
##' data.set. Then we compute a simple correlation (Pearson's by
##' default) between references predictions and the 'shuffled' 
##' one. The return score is 1-cor(pred_ref,pred_shuffled). The
##' highest the value, the more influence the variable has on the
##' model. A value of this 0 assumes no influence of that variable
##' on the model. Note that this technique does not account for
##' interactions between the variables.
##' 
##' @return a \code{list} of class "BIOMOD_variables_importances"
##' which contains:
##' 
##'  - mat: a \code{data.frame} containing variables importance
##'    scores for each permutation run.
##' 
##' @author Damien Georges
##' @seealso \code{\link[biomod2]{randomise_data}}, 
##'   \code{\link[biomod2]{full_suffling}}
##'   
##' @keywords suffle
##' @keywords random
##' @keywords importance
##' 
##' @examples
##' xx <- 
##'   data.frame( 
##'     a = sample(c(0, 1), 100, replace = TRUE),
##'     b = rnorm(100),
##'     c = 1:100
##'   )
##'   
##' mod <- glm(a ~ b + c, data = xx)
##' 
##' variables_importance(
##'   model = mod, 
##'   data = xx[, c('b', 'c')], 
##'   method = "full_rand", 
##'   nb_rand = 3
##' )
##' 
variables_importance <- 
  function(
    model, 
    data, 
    method = "full_rand", 
    nb_rand = 1, 
    ...
  ){
  out <- list()
  args <- .variables_importance.check.args(model=model, data=data, method=method, ...)

  # test prediction is computable
  ref <- try(predict(args$model, args$data))
  if(inherits(ref,"try-error")) stop("Unable to make model prediction")

  # make randomisation
  out$mat <- matrix(0,nrow=length(args$variables), ncol=nb_rand, dimnames=list(args$variables, paste('rand',1:nb_rand,sep="")))

  for(r in 1:nb_rand){
    for(v in args$variables){
      data_rand <- randomise_data(args$data,v,method)
      shuffled.pred <- predict(args$model, data_rand)
      out$mat[v,r] <- 1 - max(round(cor(x=ref, y=shuffled.pred, use="pairwise.complete.obs", method="pearson"),digits=6),0,na.rm=T)
    }
  }

  class(out) <- "BIOMOD_variables_importances"
  return(out)
}


# variables_importance argument checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
.variables_importance.check.args <- function(...){
  args <- list(...)

  # test that input data is supported
  supported_models <- c("biomod2_model", "nnet", "rpart", "fda", "gam", "glm", "lm", "gbm", "mars", "randomForest")
  if(!inherits(args$model, supported_models)) stop("Model class unsupported")

  # test method is supported
  supported_methods <- c("full_rand")
  if(! args$method %in% supported_methods ) stop("Unknown method")

  # get variables names
  if(is.null(args$variables)) args$variables <- colnames(args$data)

  return(args)
}

# data_set shuffling =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
randomise_data <- function(data,variable,method){
  if(method=='full_rand'){
    return(full_suffling(data,variable))
  }

}


# full shuffling =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
full_suffling <- function(x,id=NULL){
  if(! (is.vector(x) | is.matrix(x) | is.data.frame(x)) ) stop("x must be a 1 or 2 dimension odject")

  ## trick to ensure that the randomisation is correctly done was not the case before ##
  set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6"))*1000000)
  out <- NULL
  if(is.null(id)){
    out <- x[sample.int(length(x))]
  } else{
    out <- x
    for(idd in id){
      out[,idd] <- out[sample.int(nrow(x)),idd]
    }
  }
  return(out)
}


##### TEST #####
# setwd("~/__BIOMOD__/DevComputing/")
# x <- rbinom(n=100,size=1,prob=0.3)
# y <- rnorm(100)
# z <- rnorm(100)
# data <- as.data.frame(cbind(x,y,z))
#
# myGLM <- glm(x~y*z,family='binomial')
#
# VI <- variables_importance(myGLM, data, nb_rand=10)
#
# ##
# library(biomod2)
#
# load("GuloGulo/GuloGulo.test.models.out")
# xx <- BIOMOD_LoadModels(GuloGulo.test.models.out,models='MAXENT.Phillips')
#
# VI <- variables_importance(get(xx[1]), getModelsInputData(GuloGulo.test.models.out,'expl.var'), nb_rand=10)

#########################################
