#suppressMessages(library(foreach))
#suppressMessages(library(raster))
#suppressMessages(library(virtualspecies))

## !!!! TO DO Write the help section of the XML

#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_CNRM-ESM2-1_ssp585_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_CNRM-ESM2-1_ssp126_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_MIROC6_ssp585_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_MIROC6_ssp126_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip

#Species tables in a data collection
#Env data 1 : climate.zip --- mc85bi50.zip - mc26bi50.zip - hd85bi50.zip - hd26bi50.zip
#Env data 2 : forests.zip --- 2050_TR_LU1.zip - 2050_DP_LU1.zip

#https://pythonhosted.org/Cheetah/users_guide/errorHandling.html

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("At least 4 arguments must be supplied : \n- two input dataset files (.tabular) : metrics table and unitobs table \n- Interest variable field from metrics table \n- Response variable from unitobs table.", call. = FALSE) # if no args -> error and exit1

} else {
    source(args[1])

    ## Environmental data
    dir_env_data <- grep("env_data_", args) #Which elements from args represents environmental data ?
    n_env_data <- length(dir_env_data) #How many types of environmental data has been used ?

    zip_names <- unlist(c(args[dir_env_data + 4], strsplit(args[dir_env_data + 8], ">!<"))) #Concatenate all zip files names
    if(length(zip_names) != length(unique(zip_names))) {stop("You entered same file(s) several times in the interface or files with the same name, please check the files you entered and change the name(s) of your file(s) if necessary")} #Is there zip files with the same name or used twice in the interface ?

    if(n_env_data == 0) {stop("You need to use at least one type of environemental data to build your model")}
    n_scenario <- lapply(dir_env_data, #n_scenario : elements about future scenarios used for each type of environmental data
                                    function(x) {
                                        dir.create(args[x]) #Create a directory for each type of environmental data
                                        ## Present data
                                        unzip(args[x + 2], exdir = paste0(args[x], "/present")) #unzip present data for each type of environmental data
                                        ## Future data
                                        sfile <- strsplit(args[x + 6], ",")[[1]] #Get location of future scenarios zip files
                                        fname <- strsplit(args[x + 8], ">!<")[[1]] #Get name of future scenarios zip files
                                        sname <- gsub("\\.zip$", "", fname) #Get name of future scenarios : removing ".zip"
                                        if(length(grep("/", sname)) > 0) {sname <- gsub(".*/([^/]+)$", "\\1", sname)} #Get name of future scenarios : removing characters prior to "/" if filename is an URL 
                                        lapply(1:length(sfile),
                                               function(y) {
                                                   unzip(sfile[y], exdir = paste0(args[x], "/future/scenario_", sname[y])) #unzip each future scenarios in a dedicated directory : "env_data_x/future/scenario_y/"
                                                   if(length(list.files(paste0(args[x], "/future/scenario_", sname[y]), recursive = TRUE)) == 0) { #If unzip wasn't successful or no files were extracted
                                                       unlink(paste0(args[x], "/future/scenario_", sname[y]), recursive = TRUE) #Supress the scenario directory
                                                   }
                                               })

                                        if(length(list.files(paste0(args[x], "/present"), recursive = TRUE)) == 0 || length(list.files(paste0(args[x], "/future"), recursive = TRUE)) == 0) { #If unzip wasn't successful or no files were extracted at all for present and/or future data
                                            cat(paste0(na.omit(as.numeric(gsub("^env_data_([0-9]+)$", "\\1", args[x]))) + 1,
                                                       ": Environmental data -> unable to extract present or future archive(s), this environmental data type will be suppressed")) #Inform the user
                                            unlink(args[x], recursive = TRUE) #and supress the environmental data type directory
                                            l <- NULL
                                        } else {
                                            l <- list(args[x], #name given to the data
                                                      length(list.dirs(paste0("./", args[x], "/future"), recursive = FALSE)), #number of future scenarios used
                                                      gsub("^\\./env_data_[0-9]+/future/scenario_([^/]+)$", "\\1", list.dirs(paste0("./", args[x], "/future"), recursive = FALSE)) #names of scenarios
                                                      #length(sfile) == length(list.files(paste0("./", args[x], "/future"), recursive=T)) 
                                                      )
                                        }
                                        return(l)
                                    })

n_scenario <- n_scenario[!as.logical(lapply(n_scenario, is.null))] ## Suppress NULL values in n_scenario
if(length(n_scenario) == 0) {stop("Many of your archives were not sucessfully unzipped, please check if you used zip files and if the format 'zip' is specified in the files informations in your Galaxy history")} #If no elements in n_scenario, it means no environmental data type has both present and future data successfully unzipped

n_env_data <- lapply(n_scenario,
                     function(x) {
                         x[[1]]
                     })

###########
###########

    ## Reorganisation of files arborescence : uniformization as we use .zip files, arborescence may be different and there is several difficulties to overcome : 
    ##  - many types of raster files exists
    ##  - sometimes a raster layer is builted from a directory, sometimes from a singular file, sometimes from several files through the designation of a single file
    ##  - sometimes a raster file contains one layer (when stacked) or band (when only rasterized), sometimes it contains several layers. We consider one layer represents one environmental variable.
    ##  - sometimes it is necessary for each variables of a same environmental data type at a particular scenario to be in separated directories as several raster files are needed to build the raster and these files must have a particular name. Which makes n directories (one per variables) containing the same number of files with the same names.
    ##  - we have to check if each environmental data types has the same number of variables in the present data and in each scenarios of future data + make sure we can make they match properly

    ## First, simplify the various paths of available files : supress directories containing only directories (no files), in the end, the longest path you can have has 4 directories : One for the type of environmental data, one for present or future data, if on future data : one for each scenario and finally, when necessary, one for each variables
    files <- list.files(".", recursive = TRUE) #list available files
 
    zipf <- grep("\\.zip$", files, value = TRUE)
    if(length(zipf) > 0) {stop(".zip files have been found inside the following archives :", paste0(unique(paste0("\n- ", na.omit(as.numeric(gsub("^env_data_([0-9]+)/.*$", "\\1", zipf)))+1, ": Environmental data -> ", gsub("^env_data_[0-9]+/(present|future)/.*$", "\\1", zipf), " archive(s)"))), "\nPlease don't use nested zip files")}

    lapply(c("present","future/scenario_[^/]+"), #Wether it's present or future data
           function(x) {
               lapply(n_env_data,
                      function(y) {
                          ## We start by deleting all unnecessary directories : with a path of more directories than the maximum 4 directories prevously described
                          from_dir_long <- grep(paste0(y, "/", x, "/.+/.+/.+"), files, value = TRUE) #list files contained in more directories than only "env_data_y/x/*/" => path too long in all cases
                          to_dir_short <- gsub(paste0("(", y, "/", x, "/).+/([^/]+/[^/]+)"),"\\1\\2", from_dir_long) #write the shortened paths of those files
                          lapply(unique(gsub(paste0("(", y, "/", x, ")/.+(/[^/]+)/[^/]+"), "\\1\\2", from_dir_long)), dir.create) #create the shortened path directories
                          file.copy(from_dir_long, to_dir_short) #copy files from long to short paths
                          unlink(gsub(paste0("(", y, "/", x, "/[^/]+/).+"),"\\1", from_dir_long), recursive = TRUE) #delete unnecessary directories
                          ## Then, the last directory we kept may or may not be useful, the next steps permits to investigate its utility
                          num_dir <- length(grep(paste0("^./", y, "/", x, "/[^/]+$"), list.dirs(".", recursive = TRUE), value = TRUE)) #count the number of last directories in each paths
                          from_dir <- grep(paste0(y, "/", x, "/.+/.+"), list.files(".", recursive = TRUE), value = TRUE) #list files contained in more directories than only "env_data_y/x/"
                          to_dir <- gsub(paste0("(", y, "/", x, "/).+/([^/]+)"), "\\1\\2", from_dir) #list of paths without the last directory : where each file will eventually be copied if it is unnecessary
                          if(x == "future/scenario_[^/]+") {
                              scenario <- n_scenario[[grep(y, n_scenario)]][2] #extract the number of future scenario entered
                          }else if(x == "present") {
                              scenario <- 1
                          }
                          if(num_dir == scenario) { #if as much final directories as scenarios : last directory useless
                              file.copy(from_dir, to_dir)
                              unlink(gsub(paste0("(", y, "/", x, "/[^/]+/).+"),"\\1", from_dir), recursive = TRUE)
                          }else if(num_dir > scenario) { #if more final directories as scenarios : last directory may be also useless
                              lapply(grep(paste0("^./", y, "/", x, "/[^/]+$"), list.dirs(".", recursive = TRUE), value = TRUE),
                                     function(z) {
                                         from_dir <- list.files(z, recursive = FALSE)
                                         if(length(from_dir) == 1) { #but only if there is a single file in each directory
                                             to_dir <- gsub(paste0("(", y, "/", x, "/).+/([^/]+)"), "\\1\\2", from_dir)
                                             file.copy(from_dir, to_dir)
                                                                   }
                                     })
                          }
                      })
           })
    
    ## To overcome these issues : Try to create a raster from each directory containing files
    d <- gsub("^(.+)/[^/]+$", "\\1", list.files(".", recursive = TRUE)) #Get directories containing files
    dirs <- unique(grep(".*/.*", d, value = TRUE)) #Remove repeated directories
    #stop(paste(dirs, collapse=" "))
    rast_from_dir <- lapply(dirs, #Make a list with rasters created sucessfully from a directory
                            function(x) {
                                r <- tryCatch(raster::raster(x), error = function(e) {}) #Need to clear the  stdout/stderr : Done
                                #sink()
                                return(r)
                            })

    raster_dirs <- dirs[-grep("^NULL$", rast_from_dir)] #List of names of directories from which you can create a raster
    env_rasters <- rast_from_dir[-grep("^NULL$", rast_from_dir)] #List of rasters created from directories

    ## To overcome these issues : Try to create a raster from each environmental data files that aren't in functionning raster directories
    if(length(raster_dirs) > 0) { #If some rasters have been successfully created from a directory
        f <- grep(paste0("(", paste(raster_dirs, collapse = "|"),")"), list.files(".", recursive = TRUE)) #List those directories
        files <- list.files(".", recursive = TRUE)[-f] #Remove it from the list of files we'll try to create rasters from
    }else{
        files <- list.files(".", recursive = TRUE)
    }

    rast_from_files <- lapply(files, #Try to create a raster from each listed files
                              function(x) {
                                  r <- tryCatch(raster::raster(x), error = function(e) {}) #Need to clear the  stdout/stderr : Done
                                  return(r)
                              })

    se <- unlist(read.delim("../outputs/tool_stderr"))
    cat(paste(grep("Error in .local(.Object, ...) :", se, fixed = TRUE, invert = TRUE, value = TRUE), collapse = "\n"), "\n\n", file = "../outputs/tool_stderr", append = FALSE) ## Clear the stderr for tryCatch() errors are useless here

    raster_files <- files[-grep("^NULL$", rast_from_files)] #List of names of files from which you can create a raster 
    env_rasters <- c(env_rasters, rast_from_files[-grep("^NULL$", rast_from_files)]) #List of raster created from directories and files

    if(all(lapply(env_rasters, is.null))) {stop("The tool couldn't use any of your file as a raster")}

    names(env_rasters) <- c(raster_dirs, raster_files)
    env_rasters <- env_rasters[gtools::mixedsort(names(env_rasters))]

    ## Species data 
    ## !!!! TO DO STOP if not a tabular and if wrong colnames

    dir.create("species")
    spnb <- length(grep("--file_species", args))
    splist <- gsub(".tabular$", "", args[grep("--file_species", args) - 1])
    paste(lapply(grep("--file_species", args),
           function(x) {
               file.copy(args[x + 1], paste0("./species/", args[x - 1]))
               return(args[x - 1])
           }), "file loaded !")
    

}


###                                                 ###
### Prepare list of present and future data rasters ###
###                                                 ###

env_names <- unique(gsub("^(env_data_[0-9]+/(present|future/scenario_[^/]+))/.+$", "\\1", names(env_rasters))) ## Make a list of each type of environmental data associated with each scenario when future data

env_data_list <- lapply(env_names, #apply to each type of environmental data AND scenario if future data
                           function(x) {
                               s_env_rast <- env_rasters[grep(x, names(env_rasters))] #Select only rasters of the environmental data type AND scenario if future data x

                               r <- length(s_env_rast)
                               for(i in 1:r) { #for each raster
   
                                   if(i == 1) { #is it the first raster ? yes
                                       b <- raster::nbands(s_env_rast[[i]])
                                       if(b > 1) { #is there more than one band in the raster ? yes
                                           for(j in 1:b) { #for each band
                                                if(j == 1) { #is it the first band ? yes
                                                    env_data <- raster::stack(raster::raster(raster::filename(s_env_rast[[i]]), band = j)) #Create the stack and assuming the bands are in the natural order
                                                } else { #is it the first band ? no
                                                    env_data <- raster::addLayer(env_data, raster::raster(raster::filename(s_env_rast[[i]]), band = j))
                                                }
                                           }
                                       names(env_data) <- paste0("band_", 1:b)
                                       } else { #is there more than one band in the raster ? no
                                           env_data <- raster::stack(s_env_rast[[1]]) #Create the stack

                                       }
                                   } else { #is it the first raster ? no
                                       b <- raster::nbands(s_env_rast[[i]])

                                       if(b > 1) { #is there more than one band in the raster ? yes
                                           for(j in 1:b) { #for each band
  
                                                if(j == 1) { #is it the first band ? yes
                                                    env_data <- raster::stack(raster::raster(raster::filename(s_env_rast[[i]]), band = j)) #Create the stack assuming the bands are in the natural order

                                                } else { #is it the first band ? no
                                                    env_data <- raster::addLayer(env_data, raster::raster(raster::filename(s_env_rast[[i]]), band = j))
                                                }
                                           }
                                       names(env_data) <- paste0("band_", 1:b)
                                       } else { #is there more than one band in the raster ? no
                                           env_data <- raster::addLayer(env_data, s_env_rast[[i]])
 
                                       }
                                   }
                               }
                               names(env_data) <- paste0(x, "-", names(env_data)) ## Renaming the layers
                               return(env_data)
                           })

names(env_data_list) <- env_names

### Reduction of extent here if any ###
## !!!! TO DO Let the user choose an extent
#e <- extent(-15, 50, 40, 70)

### Reduction to the minimum extent of every layers ###

e <- lapply(env_data_list, raster::extent) ## Extract the extent of each raster

names(e) <- names(env_data_list)

comb <- as.data.frame(unique(t(apply(expand.grid(rep(list(1:length(e)), 2)), 1, sort)))) ## Find all existing 2 objects combinations
comb <- subset(comb, V1 != V2) ## Subset combination of same rasters

inte <- as.data.frame(apply(comb, 1, 
                            function(x) {
                                i <- raster::intersect(e[[x[[1]]]], e[[x[[2]]]]) ## Find whether some layers don't intersect, if so we'll have to ignore some data
                                return(c(x[[1]], x[[2]], is.null(i))) ## Return a data frame with the combinations and if they intersect : 1 if FALSE 0 if TRUE
                            }))

if(all(inte[3, ] == 1)) { ## If no rasters intersects
    env_data_list_f <- NULL
    suppr_rast <- NULL
    stop("\nNone of your rasters intersect, please check your data has intersectionnal ranges\n")

} else if(any(inte[3, ] == 1)) { ## If some raster(s) don't intersects

    suppr_rast <- c() ## Make a list for raster(s) to supress
    no_inter <- as.data.frame(inte[c(1, 2), inte[3, ] == 1]) ## Make a data frame with only combinations that doesn't intersect
  
    while(ncol(no_inter) != 0) { ## Make a loop as long as the no_inter data frame isn't empty, this loop will supress all combinations involving the raster with the least functioning intersections, if there is intersection problems remaining at the end of the first interation, it will supress all combinations involving the second raster with the least functioning intersections and so on ...
        tab <- apply(no_inter, 1, table) ## List the number of failed interactions for each raster
        max_no_inter_rast <- list2DF(lapply(tab, ## Extract names of rasters with the maximum failed interactions for each tab
                                            function(x) {
                                                return(list(names(x[max(x) == x]), as.integer(unique(x[max(x) == x]))))
                                            }))

        l <- unlist(max_no_inter_rast[1, max(unlist(max_no_inter_rast[2, ])) == max_no_inter_rast[2, ]]) ## Extract the name(s) of raster(s) with the maximum failed interactions between all rasters
        suppr_rast <- unique(unlist(c(suppr_rast, l[[1]]))) ## Add the name of the raster to supress to the suppr_rast list
        no_inter <- no_inter[, !apply(no_inter == l[[1]], 2, any)] ## Supress the combinations involving the raster to supress
    }

    env_data_list_f <- env_data_list[names(e)[-as.numeric(suppr_rast)]] ## Supress the raster(s) to supress identified in the while loop 
    cat("\nSome rasters doesn't intersect each other, they will be ignored\n Ignored raster(s) name(s): ", paste0(names(e)[as.numeric(suppr_rast)], collapse = ", "), "\n")

} else {
    env_data_list_f <- env_data_list ## No raster to supress
    suppr_rast <- NULL
    cat("\nCongrats ! All your rasters intersect each other\n")

}

if(is.null(env_data_list_f) || any(dim(env_data_list_f) == 0)) {stop("None of your rasters intersect, please check your data has intersectionnal ranges")} ## Security check

### Find min extent for the rasters ###

ext_df <- as.data.frame(do.call(rbind, lapply(names(env_data_list_f), ## Create a data frame with the extent of every raster we kept
                                              function(n) {
                                                  e <- raster::extent(env_data_list_f[[n]])
                                                  ext_l <- c(xmin = NA, xmax = NA, ymin = NA, ymax = NA)
                                                  ext_l["xmin"] <- e@xmin
                                                  ext_l["xmax"] <- e@xmax
                                                  ext_l["ymin"] <- e@ymin
                                                  ext_l["ymax"] <- e@ymax
                                                  return(ext_l)
                                              })))
rownames(ext_df) <- names(env_data_list_f)

ext_list <- lapply(colnames(ext_df), ## Create a list with the max of all xmin and ymin and the min of all xmax and ymax to get the final extent
              function(x) {
                  if(grepl("min", x)) {
                      max(as.numeric(ext_df[, x]))
                  }else if(grepl("max", x)) {
                      min(as.numeric(ext_df[, x]))
                  }
              })

ext <- raster::extent(ext_list[[1]], ext_list[[2]], ext_list[[3]], ext_list[[4]])


###              ###
### present data ###
###              ###

present_data_list <- env_data_list_f[grep("present", names(env_data_list_f))] ## Extract only present rasters

n_lay_pres <- lapply(present_data_list, raster::nlayers) ## Extract the number of layers/variables of each rasters

### Resampling of present data ###
if(raster::extent(present_data_list[[1]]) == ext) { ## If the extent of the first present raster is the same as the calculated minimum extent
    present_data_crop <- present_data_list[[1]] ## No need to crop
} else { ## If the extent of the first present raster isn't the same as the calculated minimum extent
    present_data_crop <- raster::crop(present_data_list[[1]], y = ext) ## Crop to the calculated minimum extent
}

if(length(present_data_list) > 1) { ## If there is more than one present raster
    present_data_list_rs <- lapply(present_data_list[2:length(present_data_list)], raster::resample, y = present_data_crop) ## Resample every other present rasters according to the first present raster

    present_data_list_crop <- c(present_data_crop, lapply(present_data_list_rs, ## Create a list of all rasters properly cropped and resampled
                                                          function(x) {
                                                              if(raster::extent(x) == ext) { ## If the extent of the x raster is the same as the calculated minimum extent
                                                                  data_crop <- x ## No need to crop
                                                              } else { ##  If the extent of the x raster isn't the same as the calculated minimum extent
                                                                  data_crop <- raster::crop(x, y = ext) ## Crop to the calculated minimum extent
                                                              }
                                                           return(data_crop)
                                                           }))
} else {
    present_data_list_crop <- c(present_data_crop)
}

names(present_data_list_crop) <- names(present_data_list)

### !!!!!!! Check if all future data are on same temporal series

future_data <- env_data_list_f[grep("future", names(env_data_list_f))]
future_data_list <- c()

### Check if there is the same number of layers/variables in each scenario as in present data ###

for(x in 1:length(n_env_data)) { ## For each environmental data type

     n_fut <- grep(n_env_data[x], names(future_data)) ## Extract the places of the future scenarios of the environmental data type x
     n_pres <- grep(n_env_data[x], names(n_lay_pres))  ## Extract the place of the present raster of the environmental data type x

     future_data_list_x <- future_data[n_fut]

     fut <- lapply(future_data_list_x, raster::nlayers) ## Extract the number of layers of each future scenarios of the environmental data type x
     pres <- n_lay_pres[n_pres] ## Extract the number of layers of the present raster of the environmental data type x

     pres_fut <- grep("TRUE", unlist(lapply(fut, function(y) { ## Extract the places of the future scenarios that has the same number of layers/variables as the present raster
                                                     return(y == pres)
                                                 })))

     if(length(pres_fut) == 0) { ## If none of the future scenarios has the same number of layers/variables than the present raster
         cat("\nNone of your future scenarios of environmental data type", x, "has the same number of layers/variables than the present raster, only the first layer/variable of each will be used\n")
         present_data_list_crop[[n_pres]] <- present_data_list_crop[[n_pres]][[1]] ## Subset only the first layer of the present raster
         future_data_list_x <- lapply(future_data_list_x, raster::subset, c(1)) ## Subset only the first layers of each future scenario

     } else if(length(pres_fut) != length(fut)) { ## If some of the future scenarios has the same number of layers/variables than the present raster
         cat("\nThe following future scenarios of environmental data type", x, "hasn't the same number of layers/variables than the present raster, they won't be used :\n -", paste0(gsub("^env_data_[0-9]+/future/scenario_([^/]+)$", "\\1", names(fut)[- pres_fut]), collapse = "\n - "), "\n")
         future_data_list_x <- future_data_list_x[pres_fut] ## Subset the scenarios with same number of layers

     } else if(length(pres_fut) == length(fut)) { ## If all the future scenarios has the same number of layers/variables than the present raster
         cat("\nWell done, all your future scenarios of environmental data type", x, "has the same number of layers/variables than the present raster\n")
     }
     future_data_list <- c(future_data_list, future_data_list_x)
}

cat("\n")

### Stack updated present environmental data ###

baseline <- raster::stack(virtualspecies::synchroniseNA(raster::stack(present_data_list_crop))) ## Opti ?

baseline_n <- baseline

names(baseline_n) <- gsub("^(env_data_[0-9]+).present(.?[^/]*)$", "\\1\\2", names(baseline)) ## At some point layer names ought to be the same for present and future data so deleting the "present" in the layer names here

dir.create("./data")
dir.create("./data/present")

saveRDS(baseline_n, file = paste0("./data/present/baseline_env_data_", paste0(as.numeric(gsub("^env_data_([0-9]+)/present$", "\\1", names(present_data_list_crop))) + 1, collapse = "-"), ".rds"))


###             ###
### future data ###
###             ###

### Resampling of future data ###

if(raster::extent(future_data_list[[1]]) == ext) { ## If the extent of the first future raster is the same as the calculated minimum extent
    future_data_crop <- future_data_list[[1]] ## No need to crop
} else { ## If the extent of the first future raster is different as the calculated minimum extent
    future_data_crop <- raster::crop(future_data_list[[1]], y = ext) ## Crop to the calculated minimum extent
}

if(length(future_data_list) > 1) { ## If there is more than one future raster

    future_data_list_rs <- lapply(future_data_list[2:length(future_data_list)], raster::resample, y = future_data_crop) ## Resample every other future rasters according to the first future raster, could (should?) be done for each combination with the raster of the selected scenario of the first environmental data type but WAY longer

    future_data_list_crop <- c(future_data_crop, lapply(future_data_list_rs, ## Create a list of all future rasters properly cropped and resampled
                                                          function(x) {
                                                              if(raster::extent(x) == ext) { ## If the extent of the x raster is the same as the calculated minimum extent
                                                                  data_crop <- x ## No need to crop
                                                              } else { ## If the extent of the x raster isn't the same as the calculated minimum extent
                                                                  data_crop <- raster::crop(x, y = ext) ## Crop to the calculated minimum extent
                                                              }
                                                          return(data_crop)
                                                          }))

    names(future_data_list_crop) <- names(future_data_list)

    ## Make stacks with all possible associations and let the user choose which ones are the most relevant ?

    n <- length(n_env_data) ## Number of environmental data types

    if(n == 1) { ## If only one env_data_type comb is empty so no need to stack
        future_data_stacks <- c(future_data_list_crop)

    } else {

        comb <- as.data.frame(unique(t(apply(expand.grid(rep(list(1:length(future_data_list)), n)), 1, sort)))) ## Find all existing n objects combinations

        ## Now we have to stack every combination and save it to restitute to the user so he can make a choice
        future_data_stacks <- lapply(1:nrow(comb),
                                     function(x) {
                                         s_comb <- future_data_list_crop[unlist(comb[x, ])]
                                         if(length(unique(gsub("^(env_data_[0-9]+)/future/scenario_[^/]+.*$", "\\1", names(s_comb)))) < n) { ## Combinations of rasters of the same env_data_type aren't relevant 
                                             return(NULL)
                                         } else {
                                             return(raster::stack(virtualspecies::synchroniseNA(raster::stack(s_comb))))
                                         }
                                     })

        names(future_data_stacks) <- lapply(1:nrow(comb), 
                                     function(x) {
                                         s_comb <- future_data_list_crop[unlist(comb[x, ])]
                                         if(length(unique(gsub("^(env_data_[0-9]+)/future/scenario_[^/]+.*$", "\\1", names(s_comb)))) < n) { ## Combinations of rasters of the same env_data_type aren't relevant 
                                             return(x)
                                         } else {
                                             return(paste0(gsub("^env_data_[0-9]+/future/scenario_([^/]+.*)$", "\\1", names(s_comb)), collapse = "-"))
                                         }
                                     })
    }

} else {
    future_data_stacks <- c(future_data_crop)
}

future_data_stacks <- future_data_stacks[!as.logical(lapply(future_data_stacks, is.null))] ## Suppress NULL values in future_data_stacks

dir.create("./data/future")

lapply(names(future_data_stacks), function(x) {
                               names(future_data_stacks[[x]]) <- names(baseline)
                               saveRDS(future_data_stacks[x], file = paste0("./data/future/scenarios_combination_", x, ".rds"))
                           })

## !!!! TO DO Verify if saved data .rds is working on R and QGIS
## !!!! TO DO Diagnostic plots
#dir.create("graphiques")
#png("./graphiques/diff_bio1.png", width = 600, height = 600)
#plot(futureData[["bio1"]] - baseline[["bio1"]])





####                          ####
#### Section 2 : Species data ####
####                          ####

## !!!! TO DO Translate and comment more on this part
## !!!! TO DO Diagnostic plots

#baseline <- stack("./data/baseline")
lapply(splist, 
       function(sp) {
  sp.occurrence <- read.table(paste0("./species/", sp, ".tabular"), 
                              sep = "\t", h = T)
  #plot(baseline[[1]])
  #points(sp.occurrence[, c("x", "y")], cex = .5, 
   #      pch = c(1, 16)[sp.occurrence$Observed + 1]) ## Visualisation of the data to make at the end
  
  # 1. Rasterisation des occurrences
  # On transforme le champ "Observed" (compose de 1 & 0) en somme par pixel
  sp.env <- raster::rasterize(x = sp.occurrence[, c("x", "y")], 
                              y = baseline,
                              field = sp.occurrence[, "Observed"],
                              fun = function(x, ...) sum(x, ...)) # Bien verifier le resultat et ajuster la fonction
  
  #plot(sp.env, col = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000")) ## SAME make it at the end
  sp.env[sp.env > 1] <- 1 ## Transform abundance into presence-absence
  # On ajoute le nom de l'espece au raster
  names(sp.env) <- sp

  # On ajoute les donnees environnementales aux presences rasterisees
  #stop(sp.env))
  sp.env <- raster::stack(sp.env,
                          baseline)

  # 2. Suppression des presences pour lesquelles on n'a pas de donnees environnementales
  # Recuperation des coordonnees de toutes les cellules
  coorXY <- raster::xyFromCell(baseline, 1:raster::ncell(baseline))
  # Transformation du raster en data.frame pour tester les NAs
  sp.env.df <- raster::getValues(sp.env)
  
  # Introduction volontaire d'erreurs pour l'exemple
  # sp.env.df[which(is.na(sp.env.df[, "bio1"])), 1] <- 1


  ## Pq bio1 juste ??? Parce que les NA sont synchro ?


  if(any(is.na(sp.env.df[, names(baseline)[1]]) & !is.na(sp.env.df[, sp])))
  {
    cat("Some points are in pixels without environmental values\n")
  }

  # On supprime les cellules pour lesquelles on n'a pas de donnees environnementales
  coorXY <- coorXY[-which(is.na(sp.env.df[, names(baseline)[1]])), ]
  sp.env.df <- sp.env.df[-which(is.na(sp.env.df[, names(baseline)[1]])), ]

  # Nombre de cellules resultantes :
  cat(sp, "\nNumber of pixels of presence:",
      "\n - Initial: ", length(which(sp.occurrence$Observed == 1)),
      "\n - After rasterisation: ", length(which(sp.env.df[, 1] == 1)), "\n")
  cat(sp, "\nNumber of pixels of absence:",
      "\n - Initial: ", length(which(sp.occurrence$Observed == 0)),
      "\n - After rasterisation: ", length(which(sp.env.df[, 1] == 0)), "\n\n")

  # 3. Recuperation des occurrences rasterisees et ecriture sur le disque
  P.points <- data.frame(
    # D'abord on recupere les coordonnees XY qui correspondent a nos cellules de presences/absences
    coorXY[which(!is.na(sp.env.df[, sp])), ],
    # Ensuite, on recupere la colonne qui indique presence/absence pour chaque cellule
    Observed = sp.env.df[which(!is.na(sp.env.df[, 1])), sp]) # On recupere les occurrences ici
  saveRDS(P.points, file = paste0("./data/occurrences_", sp))
})

