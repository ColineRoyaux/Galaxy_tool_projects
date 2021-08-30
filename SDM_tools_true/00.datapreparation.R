suppressMessages(library(raster))
#suppressMessages(library(virtualspecies))
suppressMessages(library(rgdal))

#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_CNRM-ESM2-1_ssp585_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_CNRM-ESM2-1_ssp126_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_MIROC6_ssp585_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut/10m/wc2.1_10m_bioc_MIROC6_ssp126_2021-2040.zip
#https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip
#https://pythonhosted.org/Cheetah/users_guide/errorHandling.html

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("At least 4 arguments must be supplied : \n- two input dataset files (.tabular) : metrics table and unitobs table \n- Interest variable field from metrics table \n- Response variable from unitobs table.", call. = FALSE) # if no args -> error and exit1

} else {
    source(args[1])

    ## Environmental data
    dir_env_data <- grep("env_data_", args) #Which elements from args represents environmental data ?
    n_env_data <- length(dir_env_data) #How many types of environmental data has been used ?
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
    ##  - sometimes a raster file contains one layer (or "band"), sometimes it contains several layers. We consider one layer represents one environmental variable.
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
                                r <- tryCatch(raster(x), error = function(e){}) #Need to clear the  stdout/stderr
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
                                  r <- tryCatch(raster(x), error = function(e) {}) #Need to clear the  stdout/stderr
                                  return(r)
                                          })
    raster_files <- files[-grep("^NULL$", rast_from_files)] #List of names of files from which you can create a raster 
    env_rasters <- c(env_rasters, rast_from_files[-grep("^NULL$", rast_from_files)]) #List of raster created from directories and files
    #stop(length(raster_files))
    #stop(length(grep("env_data_0", raster_files, value=TRUE)))
    #stop(paste(grep("env_data_0", raster_files, value=TRUE), collapse = " "))

    ## Species data
    dir.create("species")
    spnb <- length(grep("--file_species", args))
    splist <- gsub(".tabular$", "", args[grep("--file_species", args)-1])
    lapply(grep("--file_species", args),
           function(x) {
               file.copy(args[x+1], paste0("./species/", args[x-1]))
               return(paste(args[x-1], "file loaded !"))
                       })

}

#### Section 1 : donnees environnementales ####
# Verifier la compatibilite des donnees ! Georeferencement, etc.

#env_data_trt <- function(dir) {
#setwd(dir)


#}

# 1. Donnees baseline (worldclim 1950-2000)
# 1.1 Donnees climatiques worldclim

climate.vars <- list.files("./env_data_0/present/", full.names=FALSE, recursive=FALSE)
climate.vars <- climate.vars[-length(climate.vars)]

for (j in climate.vars)
{
  if(j == climate.vars[1]) # Pour le georeferencement
  {
    envData <- stack(raster(paste0("./env_data_0/present/", j)))
  } else
  {
    envData <- addLayer(envData, 
                        raster(paste0("./env_data_0/present/", j)))
  }
}

names(envData) <- paste0("bio", 1:19) # Pour avoir les mêmes noms que les fichiers futurs
#names(envData) <- gsub("wc2.1_10m_bio_([0-9]+).tif","bio\\1",climate.vars) # Pour avoir les mêmes noms que les fichiers futurs

# Reduction a l'ouest-palearctique
e <- extent(-15, 50, 40, 70)
envData <- crop(envData, e)

# 1.2 Autres donneesbaseline 
forests <- raster("./env_data_1/present/2000_LU1.grd")

names(forests) <- "forests"

# Reechantillonnage : aggregate, disaggregate (attention !), resample
forests.0.167 <- resample(x = forests, 
                          y = envData)
forests.0.167 <- crop(forests.0.167, 
                      envData)
envData <- stack(envData,
                 forests.0.167)

# Synchronisation des NAs : explication
val <- getValues(envData)
which(is.na(val[, "bio1"]) & !is.na(val[, "forests"]))
length(which(is.na(val[, "forests"]) & !is.na(val[, "bio1"])))
dummy.raster <- envData[[1]]
dummy.raster[] <- NA
dummy.raster[which(is.na(val[, "bio1"]) & !is.na(val[, "forests"]))] <- 1
dummy.raster[which(is.na(val[, "forests"]) & !is.na(val[, "bio1"]))] <- 2
plot(dummy.raster, col = c("red", "blue"))

envData <- synchroniseNA(envData)

# Ecriture des donnees
dir.create("data")
writeRaster(envData, "./data/baseline", overwrite = T) # Overwrite?

# 2. Donnees futures
#future.datasets <- list.files("./env_data_0/future/", full.names = TRUE, recursive = FALSE)
future.datasets <- c("hd26bi50", # HadGEM2-AO, RCP 2.6
                     "hd85bi50", # HadGEM2-AO, RCP 8.5
                     "mc26bi50", # MIROC5, RCP 2.6
                     "mc85bi50") # MIROC5, RCP 8.5

future.names <- c("CNRM-ESM2_RCP2.6_2040",
                  "CNRM-ESM2_RCP8.5_2040",
                  "MIROC6_RCP2.6_2040",
                  "MIROC6_RCP8.5_2040")

baseline <- stack("./data/baseline")

# Important : les noms de variables doivent être identiques entre baseline et futur
i <- 0
for(gcm in future.datasets)
{
  i <- i + 1
  cat(future.names[i], "\n")
  
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
      stop("lol")
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
