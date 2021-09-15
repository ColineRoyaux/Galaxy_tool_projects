#suppressMessages(library(raster))
#suppressMessages(library(virtualspecies))
#suppressMessages(library(rgdal))

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
    if(n_env_data == 0) {stop("You need to use at least one type of environemental data to build your model")}
    n_scenario <- lapply(dir_env_data, #n_scenario : elements about future scenarios used for each type of environmental data
                         function(x) {
                             dir.create(args[x]) #Create a directory for each type of environmental data
                             ## Present data
                             unzip(args[x+2], exdir = paste0(args[x], "/present")) #unzip present data for each type of environmental data
                             ## Future data
                             sfile <- strsplit(args[x+4], ",")[[1]] #Get location of future scenarios zip files
                             fname <- strsplit(args[x+6], ">!<")[[1]] #Get name of future scenarios zip files
                             sname <- gsub("\\.zip$", "", fname) #Get name of future scenarios : removing ".zip"
                             if(length(grep("/", sname)) > 0) {sname <- gsub(".*/([^/]+)$", "\\1", sname)} #Get name of future scenarios : removing characters prior to "/" if filename is an URL 
                             #if(args[x] == "env_data_1"){(paste(sname, collapse=" "))}
                             lapply(1:length(sfile),
                                    function(y) {
                                        unzip(sfile[y], exdir = paste0(args[x], "/future/scenario_", sname[y])) #unzip each future scenarios in a dedicated directory : "env_data_x/future/scenario_y/"
                                    })
                             l <- list(args[x], #name given to the data
                                       length(sfile), #number of future scenarios used
                                       sname #names of scenarios
                                       #length(sfile) == length(list.files(paste0("./", args[x], "/future"), recursive=T)) 
                                       )
                             return(l)
                         })

###########
########### ADD A LINE TO UNZIP ZIP FILES INSIDE THE ZIP => HAPPENS SOMETIMES TO HAVE ZIP ARCHIVE(S) INSIDE A ZIP ARCHIVE 
###########

    ## Reorganisation of files arborescence : uniformization as we use .zip files, arborescence may be different and there is several difficulties to overcome : 
    ##  - many types of raster files exists
    ##  - sometimes a raster layer is builted from a directory, sometimes from a singular file, sometimes from several files through the designation of a single file
    ##  - sometimes a raster file contains one layer (when stacked) or band (when only rasterized), sometimes it contains several layers. We consider one layer represents one environmental variable.
    ##  - sometimes it is necessary for each variables of a same environmental data type at a particular scenario to be in separated directories as several raster files are needed to build the raster and these files must have a particular name. Which makes n directories (one per variables) containing the same number of files with the same names.
    ##  - we have to check if each environmental data types has the same number of variables in the present data and in each scenarios of future data + make sure we can make they match properly

    ## First, simplify the various paths of available files : supress directories containing only directories (no files), in the end, the longest path you can have has 4 directories : One for the type of environmental data, one for present or future data, if on future data : one for each scenario and finally, when necessary, one for each variables
    files <- list.files(".", recursive = TRUE) #list available files
    lapply(c("present","future/scenario_[^/]+"), #Wether it's present or future data
           function(x) {
               lapply(args[dir_env_data],
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
                          to_dir <- gsub(paste0("(", y, "/", x, "/).+/([^/]+)"),"\\1\\2", from_dir) #list of paths without the last directory : where each file will eventually be copied if it is unnecessary
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
                                             to_dir <- gsub(paste0("(", y, "/", x, "/).+/([^/]+)"),"\\1\\2", from_dir)
                                             file.copy(from_dir, to_dir)
                                                                   }
                                     })
                          }
                      })
           })

    #stop(paste(grep("present", list.files(".", recursive = TRUE), value=TRUE), collapse = " "))
   
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

    #stop(class(rast_from_dir[[5]]))
    raster_dirs <- dirs[-grep("^NULL$", rast_from_dir)] #List of names of directories from which you can create a raster
    env_rasters <- rast_from_dir[-grep("^NULL$", rast_from_dir)] #List of rasters created from directories
    #stop(length(raster_dirs))
    #stop(paste(raster_dirs, collapse = " "))

    ## To overcome these issues : Try to create a raster from each environmental data files that aren't in functionning raster directories
    if(length(raster_dirs) > 0) { #If some rasters have been successfully created from a directory
        f <- grep(paste0("(", paste(raster_dirs, collapse = "|"),")"), list.files(".", recursive = TRUE)) #List those directories
        files <- list.files(".", recursive = TRUE)[-f] #Remove it from the list of files we'll try to create rasters from
    }else{
        files <- list.files(".", recursive = TRUE)
    }
    #stop(length(files))
    #stop(paste(files, collapse = " "))

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
    #stop(grep("future", names(env_rasters), value = TRUE))
    #stop(length(raster_files))
    #stop(length(grep("env_data_0", raster_files, value=TRUE)))
    #stop(paste(grep("env_data_0", raster_files, value=TRUE), collapse = " "))

    ## Species data
    dir.create("species")
    spnb <- length(grep("--file_species", args))
    splist <- gsub(".tabular$", "", args[grep("--file_species", args)-1])
    paste(lapply(grep("--file_species", args),
           function(x) {
               file.copy(args[x+1], paste0("./species/", args[x-1]))
               return(args[x-1])
           }), "file loaded !")
    

}

###                                                 ###
### Prepare list of present and future data rasters ###
###                                                 ###

env_names <- unique(gsub("^(env_data_[0-9]+/(present|future/scenario[^/]+))/.+$", "\\1", names(env_rasters))) ## Made a list of each type of environmental data associated with each scenario when future data

env_data_list <- lapply(env_names, #apply to each type of environmental data AND scenario if future data
                           function(x) {
                               s_env_rast <- env_rasters[grep(x, names(env_rasters))] #Select only rasters of the environmental data type AND scenario if future data x

                               r <- length(s_env_rast)
                               for(i in 1:r) { #for each raster
   
                                   if(i == 1) { #is it the first raster ? yes
                                       b <- raster::nbands(s_env_rast[[i]])

                                       if(b > 1) { #is there more than one band in the raster ? yes !!!!!!! CONDITION TO TEST !!!!!!
                                           stop("data with several bands, test on !!!!")
                                           for(j in 1:b) { #for each band

                                                if(j == 1) { #is it the first band ? yes
                                                    env_data <- raster::stack(raster::raster(names(s_env_rast[[i]]), band = j)) #Create the stack

                                                } else { #is it the first band ? no
                                                    env_data <- raster::addLayer(env_data, raster::raster(names(s_env_rast[[i]]), band = j))
                                                }
                                           }
  
                                       } else { #is there more than one band in the raster ? no
                                           env_data <- raster::stack(s_env_rast[[1]]) #Create the stack

                                       }
                                   } else { #is it the first raster ? no
                                       b <- raster::nbands(s_env_rast[[i]])

                                       if(b > 1) { #is there more than one band in the raster ? yes !!!!!!! CONDITION TO TEST !!!!!!
                                           stop("data with several bands, test on !!!!")
                                           for(j in 1:b) { #for each band
  
                                                if(j == 1) { #is it the first band ? yes
                                                    env_data <- raster::stack(raster::raster(names(s_env_rast[[i]]), band = j)) #Create the stack

                                                } else { #is it the first band ? no
                                                    env_data <- raster::addLayer(env_data, raster::raster(names(s_env_rast[[i]]), band = j))
                                                }
                                           }

                                       } else { #is there more than one band in the raster ? no
                                           env_data <- raster::addLayer(env_data, s_env_rast[[i]])
 
                                       }
                                   }
                               }
                               names(env_data) <- paste0(x, "-", names(env_data)) ## Renaming the layers, present/future necessary ?
                               return(env_data)
                           })

names(env_data_list) <- env_names

stop(names(env_data_list))

### Reduction of extent here if any ###
#e <- extent(-15, 50, 40, 70)

### Reduction to the minimum extent of every layers ###

#apply(ext_df, 1, function(x) {stop(x)})
#example data :
#env_data_list <- list(E = raster::extent(10, 14 , -3, 1), A = raster::extent(-16, 0, -11, -5), C = raster::extent(-9, 0, -12, 10), D = raster::extent(-7, 17, -9, 11), B = raster::extent(-17, 8, -12, 8), G = raster::extent(-19, 19, -13, 13))

e <- lapply(env_data_list, raster::extent)

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
    stop("None of your rasters intersect, please check your data has intersectionnal ranges")

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
    cat("Some rasters doesn't intersect each other, they will be ignored\n Ignored raster(s) name(s): ", paste0(names(e)[as.numeric(suppr_rast)], collapse = ", "))

} else {
    env_data_list_f <- env_data_list ## No raster to supress
    suppr_rast <- NULL
    cat("Congrats ! All your rasters intersect each other")

}











###                                     ###
### Baseline data "present" directories ###
###                                     ###



### Reduction of extent here if any ###
#e <- extent(-15, 50, 40, 70)

### Reduction to the minimum extent of every layers ###

#apply(ext_df, 1, function(x) {stop(x)})
#example data :
#env_data_list <- list(E = raster::extent(10, 14 , -3, 1), A = raster::extent(-16, 0, -11, -5), C = raster::extent(-9, 0, -12, 10), D = raster::extent(-7, 17, -9, 11), B = raster::extent(-17, 8, -12, 8), G = raster::extent(-19, 19, -13, 13))

e <- lapply(env_data_list, raster::extent)

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
    stop("None of your rasters intersect, please check your data has intersectionnal ranges")

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
    cat("Some rasters doesn't intersect each other, they will be ignored\n Ignored raster(s) name(s): ", paste0(names(e)[as.numeric(suppr_rast)], collapse = ", "))

} else {
    env_data_list_f <- env_data_list ## No raster to supress
    suppr_rast <- NULL
    cat("Congrats ! All your rasters intersect each other")

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

### Resampling of data ###

env_data_crop <- raster::crop(env_data_list_f[[1]], y = ext)
#env_data_list_f[[1]] <- raster::brick(raster::crop(env_data_list_f[[1]], y = ext), values = TRUE, nl = raster::nlayers(env_data_list_f[[1]]))

env_data_list_f_rs <- lapply(env_data_list_f[2:length(env_data_list_f)], function(n) {
  return(raster::resample(x = n, y = env_data_list_f[[1]]))
})

env_data_list_crop <- c(env_data_crop, lapply(env_data_list_f_rs, raster::crop, y = ext))

### Stack present environmental data ###

baseline <- raster::stack(virtualspecies::synchroniseNA(raster::stack(env_data_list_crop))) ## Opti ?



###             ###
### future data ###
###             ###

future_datasets <- lapply(n_scenario, 
                          function(x) {
                              return(unlist(x[3]))
                          })
names(future_datasets) <- unlist(lapply(n_scenario, 
                                 function(x) {
                                     return(x[1])
                                 }))
#fut_env_scenario <- unlist(lapply(args[dir_env_data], #apply to each type of environmental data
 #                          function(x) {
  #                             lapply(future_datasets[x], #apply to each scenario
   #                                   function(y) {
    #                                      return(paste0(x, "/future/scenario_", y))
     #                                 })
      #                     })) ## Made a list of each type of environmental data associated with each scenario


future_data_list_stack <- lapply(args[dir_env_data],
                                 function(x) {
                                     return(raster::stack(future_data_list[grep(x, names(future_data_list))]))
                                 })

stop(names(future_data_list_stack[[2]]))

### Resampling of data ### ## Make a function ?? Approx the same as present data 

future_data_crop <- raster::crop(future_data_list_stack[[1]], y = ext)

future_data_list_rs <- lapply(future_data_list_stack[2:length(env_data_list_f)], function(n) {
  return(raster::resample(x = n, y = env_data_list_f[[1]]))
})

future_data_list_crop <- c(env_data_crop, lapply(env_data_list_f_rs, raster::crop, y = ext))

### Stack future environmental data ###

baseline <- raster::stack(virtualspecies::synchroniseNA(raster::stack(env_data_list_crop))) ## Opti ?

#i <- 0
#for(gcm in future_datasets)
#{
  #i <- i + 1
  #cat(future.names[i], "\n")
  
  # Donnees climatiques
  
  #for (j in 1:nbands(raster(paste0("./env_data_0/future/", gcm))))
  #{
    #if(j == 1)
    #{
    #  futureData <- stack(raster(paste0(gcm, band = j))
#    } else {
 #     futureData <- addLayer(futureData, raster(paste0("./env_data_0/future/", gcm), band = j))
  #  }
  #}
  #stop(nlayers(futureData))

  #climate.vars <- paste0(gcm, 1:19)
  #stop(paste(list.files("./env_data_1/future/", recursive = TRUE), collapse = " "))
 # for (j in climate.vars)
#  {
    #if(j == climate.vars[1])
   # {
     # stop(paste0("./env_data_0/future/scenario_", gcm, "/", j, ".tif"))
  #    futureData <- stack(raster(paste0("./env_data_0/future/scenario_", gcm, "/", j, ".tif")))
 #   } else
#    {
    #  futureData <- addLayer(futureData, raster(paste0("./env_data_0/future/scenario_", gcm, "/", j, ".tif")))
   # }
  #}
  #futureData <- crop(futureData, baseline)

 # names(futureData) <- paste0("bio", 1:19)

  # Autres donnees ## ????????????????? Assigner avec l'interface ??????????????
  if(length(grep("26", gcm)))
  {
    LU.scenario <- "DP"
  } else if(length(grep("85", gcm)))
  {
    LU.scenario <- "TR"
  }
  #forests <- raster(paste0("./env_data_1/future/scenario_2050_", LU.scenario, "_LU1/", "2050_", LU.scenario, "_LU1.grd")) # Programmer annee si plusieurs periodes futures
  #names(forests) <- "forests"
    #stop("lol")
  forests.167 <- resample(forests, 
                          futureData)
  forests.167 <- crop(forests.167, futureData)
  
  futureData <- stack(futureData,
                      forests.167)

  # Synchronisation des stacks futurs + synchronisation par rapport au baseline
  tmp <- baseline[["bio1"]]
  names(tmp) <- "tmp"
  futureData <- addLayer(tmp, futureData)
  futureData <- synchroniseNA(futureData)
  futureData <- dropLayer(futureData, 1)
  
  writeRaster(futureData, paste("./data/", future.names[i], ".grd", sep = ""), overwrite = T) 
  saveRDS(futureData, 
       file = paste0("./data/", future.names[i], ".rds"))
  cat(paste0(Sys.time(), " - ", future.names[i], "  done\n"))
  
}



#Infos sur les scénarios futurs des différentes données env dans n_scenario
#Trouver comment extraire les noms de variables : a partir des noms de dossier si plusieurs fichiers par variable et a partir des noms de fichier sinon MAIS soucis parce que parfois c'est en un seul fichier on a plusieurs couches qui correspondent aux différentes variables

   # mult_var <- lapply(names(env_rasters),
    #            function(x) {
     #               if(grepl("present", x) && length(unlist(gregexpr("/", x))) > 1) {
      #                  return
       #             }else if(grepl("future", x)) {
        #                return(length(unlist(gregexpr("/", x))) > 2)
         #           }
          #                  })
    #stop(mult_var)

#### Section 1 : donnees environnementales ####
# Verifier la compatibilite des donnees ! Georeferencement, etc.

#env_data_trt <- function(dir) {
#setwd(dir)


#}

# 1. Donnees baseline (worldclim 1950-2000)
# 1.1 Donnees climatiques worldclim

#climate.vars <- list.files("./env_data_0/present/", full.names=FALSE, recursive=FALSE)

#climate.vars <- climate.vars[-length(climate.vars)]

#for (j in climate.vars)
#{
 # if(j == climate.vars[1]) # Pour le georeferencement
  #{
   # envData <- stack(raster(paste0("./env_data_0/present/", j)))
#  } else {
 #   envData <- addLayer(envData, 
  #                      raster(paste0("./env_data_0/present/", j)))
#  }
#}

#names(envData) <- paste0("bio", 1:19) # Pour avoir les mêmes noms que les fichiers futurs
#names(envData) <- gsub("wc2.1_10m_bio_([0-9]+).tif","bio\\1",climate.vars) # Pour avoir les mêmes noms que les fichiers futurs

# Reduction a l'ouest-palearctique
#e <- extent(-15, 50, 40, 70)
#envData <- crop(envData, e)

# 1.2 Autres donneesbaseline 
#forests <- raster("./env_data_1/present/2000_LU1.grd")

#names(forests) <- "forests"

# Reechantillonnage : aggregate, disaggregate (attention !), resample
#forests.0.167 <- resample(x = forests, 
     #                     y = envData)
#forests.0.167 <- crop(forests.0.167, 
      #                envData)
#envData <- stack(envData,
       #          forests.0.167)

## Synchronisation des NAs : explication
#val <- getValues(envData)
#which(is.na(val[, "bio1"]) & !is.na(val[, "forests"]))
#length(which(is.na(val[, "forests"]) & !is.na(val[, "bio1"])))
#dummy.raster <- envData[[1]]
#dummy.raster[] <- NA
#dummy.raster[which(is.na(val[, "bio1"]) & !is.na(val[, "forests"]))] <- 1
#dummy.raster[which(is.na(val[, "forests"]) & !is.na(val[, "bio1"]))] <- 2
#plot(dummy.raster, col = c("red", "blue"))

#envData <- synchroniseNA(envData)

# Ecriture des donnees
#dir.create("data")
#writeRaster(envData, "./data/baseline", overwrite = T) # Overwrite?

# 2. Donnees futures
#future.datasets <- list.files("./env_data_0/future/", full.names = TRUE, recursive = FALSE)
#future.datasets <- c("hd26bi50", # HadGEM2-AO, RCP 2.6
 #                    "hd85bi50", # HadGEM2-AO, RCP 8.5
  #                   "mc26bi50", # MIROC5, RCP 2.6
   #                  "mc85bi50") # MIROC5, RCP 8.5

#future.names <- c("CNRM-ESM2_RCP2.6_2040",
 #                 "CNRM-ESM2_RCP8.5_2040",
  #                "MIROC6_RCP2.6_2040",
   #               "MIROC6_RCP8.5_2040")

#baseline <- stack("./data/baseline")

# Important : les noms de variables doivent être identiques entre baseline et futur
i <- 0
for(gcm in future.datasets)
{
  i <- i + 1
  #cat(future.names[i], "\n")
  
  # Donnees climatiques
  
  #for (j in 1:nbands(raster(paste0("./env_data_0/future/", gcm))))
  #{
    #if(j == 1)
    #{
    #  futureData <- stack(raster(paste0(gcm, band = j))
#    } else {
 #     futureData <- addLayer(futureData, raster(paste0("./env_data_0/future/", gcm), band = j))
  #  }
  #}
  #stop(nlayers(futureData))

  climate.vars <- paste0(gcm, 1:19)
  #stop(paste(list.files("./env_data_1/future/", recursive = TRUE), collapse = " "))
  for (j in climate.vars)
  {
    if(j == climate.vars[1])
    {
     # stop(paste0("./env_data_0/future/scenario_", gcm, "/", j, ".tif"))
      futureData <- stack(raster(paste0("./env_data_0/future/scenario_", gcm, "/", j, ".tif")))
    } else
    {
      futureData <- addLayer(futureData, raster(paste0("./env_data_0/future/scenario_", gcm, "/", j, ".tif")))
    }
  }
  futureData <- crop(futureData, baseline)

  names(futureData) <- paste0("bio", 1:19)

  # Autres donnees
  if(length(grep("26", gcm)))
  {
    LU.scenario <- "DP"
  } else if(length(grep("85", gcm)))
  {
    LU.scenario <- "TR"
  }
  forests <- raster(paste0("./env_data_1/future/scenario_2050_", LU.scenario, "_LU1/", "2050_", LU.scenario, "_LU1.grd")) # Programmer annee si plusieurs periodes futures
  names(forests) <- "forests"
    #stop("lol")
  forests.167 <- resample(forests, 
                          futureData)
  forests.167 <- crop(forests.167, futureData)
  
  futureData <- stack(futureData,
                      forests.167)

  # Synchronisation des stacks futurs + synchronisation par rapport au baseline
  tmp <- baseline[["bio1"]]
  names(tmp) <- "tmp"
  futureData <- addLayer(tmp, futureData)
  futureData <- synchroniseNA(futureData)
  futureData <- dropLayer(futureData, 1)
  
  writeRaster(futureData, paste("./data/", future.names[i], ".grd", sep = ""), overwrite = T) 
  saveRDS(futureData, 
       file = paste0("./data/", future.names[i], ".rds"))
  cat(paste0(Sys.time(), " - ", future.names[i], "  done\n"))
  
}

# Synchronisation du stack baseline par rapport aux futurs
tmp <- futureData[["bio1"]]
names(tmp) <- "tmp"

baseline <- addLayer(tmp,
                     baseline)
baseline <- synchroniseNA(baseline)
baseline <- dropLayer(baseline, 1)
baseline <- stack(baseline) # to remove brick status
writeRaster(baseline, 
            "./data/baseline", overwrite = T)
saveRDS(baseline, 
     file = "./data/baseline.rds")

png("./data/diff_bio1.png", width = 600, height = 600)
plot(futureData[["bio1"]] - baseline[["bio1"]])
dev.off()



#### Section 2 : donnees espece ####

baseline <- stack("./data/baseline")
for (sp in splist)
{
  sp.occurrence <- read.table(paste0("./species/", sp, ".tabular"), 
                              sep = "\t", h = T)
  plot(baseline[[1]])
  points(sp.occurrence[, c("x", "y")], cex = .5, 
         pch = c(1, 16)[sp.occurrence$Observed + 1])
  
  # 1. Rasterisation des occurrences
  # On transforme le champ "Observed" (compose de 1 & 0) en somme par pixel
  sp.env <- rasterize(x = sp.occurrence[, c("x", "y")], 
                      y = baseline,
                      field = sp.occurrence[, "Observed"],
                      fun = function(x, ...) sum(x, ...)) # Bien verifier le resultat et ajuster la fonction
  
  plot(sp.env, col = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"))
  sp.env[sp.env > 1] <- 1
  # On ajoute le nom de l'espece au raster
  names(sp.env) <- sp

  # On ajoute les donnees environnementales aux presences rasterisees
  sp.env <- stack(sp.env,
                  baseline)

  # 2. Suppression des presences pour lesquelles on n'a pas de donnees environnementales
  # Recuperation des coordonnees de toutes les cellules
  coorXY <- xyFromCell(baseline, 1:ncell(baseline))
  # Transformation du raster en data.frame pour tester les NAs
  sp.env.df <- getValues(sp.env)
  
  # Introduction volontaire d'erreurs pour l'exemple
  # sp.env.df[which(is.na(sp.env.df[, "bio1"])), 1] <- 1
  
  if(any(is.na(sp.env.df[, "bio1"]) & !is.na(sp.env.df[, sp])))
  {
    cat("Some points are in pixels without environmental values\n")
  }

  # On supprime les cellules pour lesquelles on n'a pas de donnees environnementales
  coorXY <- coorXY[-which(is.na(sp.env.df[, "bio1"])), ]
  sp.env.df <- sp.env.df[-which(is.na(sp.env.df[, "bio1"])), ]
  
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
}
