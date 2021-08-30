library(biomod2)
library(RColorBrewer)
library(raster)

setwd("/home/pndb-cr/2-Galaxy/SDM/Test_input")

# rasterOptions(tmpdir = "d:/raster_temp")
# Attention à charger le stack avec le répertoire complet si
# vous utilisez stack()
baseline <- readRDS("./data/baseline.rds")

sp.list <- read.csv("./data_cours/species_list.csv", sep = ";")

listStacks <- c("baseline",
                "HadGEM2AO_RCP2.6_2050",
                "HadGEM2AO_RCP8.5_2050",
                "MIROC5_RCP2.6_2050",
                "MIROC5_RCP8.5_2050")


# Une fois que la sélection de variable sera effectuée, on utilisera cette commande
# load("./data/selected.variables.RData")
# En attendant, on créer une liste avec les mêmes variables pour tout le monde:
sel.vars <- lapply(as.character(sp.list$sp),
                   function (x) c("bio5", "bio6", "forests"))
names(sel.vars) <- as.character(sp.list$sp)


# On stocke le répertoire de travail initial
initial.wd <- getwd()

for (i in 1:nrow(sp.list))
{
  # Pour tester la boucle:
  # i <- 3
  
  # Nom de l'espèce
  sp <- as.character(sp.list$sp[i])

  # Points d'occurrence de l'espèce
  P.points <- readRDS(paste0("./data/occurrences_", sp))  
  
  # Variables sélectionnes pour l'espèce
  cur.vars <- sel.vars[[sp]]
  
  # 1. Biomod2 doit-il générer des pseudoabsences ?
  if(sp.list$pa.generation[i] == "biomod")
  {
    nb.PA <- nrow(P.points) # Générer autant de pseudoabsences
    # que de présences
    runs.PA <- 3 # Nombre de répétitions de pseudo-absences
  } else
  {
    runs.PA <- 0
    nb.PA <- 0
  }
  
  # 2. Les données environnementales viennent d'un raster ou d'une matrice ?
  # Dans votre cas actuel tout est stocké en objet spatial (i.e. rasters)
  if(sp.list$env.data.type[i] == "Spatial")
  { 
    calib.env.data <- baseline[[cur.vars]]
  } else
  # Cas où les données environnementales sont stockées sous forme de matrice
  {
    calib.env.data <- P.points[, cur.vars]
  }
  
  # 3. Formatage des occurrences pour biomod
  # On crée un objet avec les coordonnées XY
  coorxy <- P.points[, c("x", "y")]
  # Et un objet avec l'information présence / absence
  P.points <- P.points[, "Observed"]
  
  
  # 4. Préparation d'un sous-dossier pour notre espèce
  if(!exists(paste0(initial.wd, "/models/", sp)))
  {
    dir.create(paste0(initial.wd, "/models/", sp), recursive = T)
  }
  setwd(paste0(initial.wd, "/models/"))
  
  
  # 5. Fonction d'initialisation de biomod
  run.data <- BIOMOD_FormatingData(resp.var = P.points, # Variable réponse (occurrences)
                                   expl.var = calib.env.data, # Variables explicatives
                                   resp.name = sp, # Nom de l'espèce
                                   resp.xy = coorxy, # Coordonnées d'occurrence
                                   PA.nb.rep = runs.PA, # Nombre de runs de pseudoabsence
                                   PA.nb.absences = nb.PA, # Nombre de pseudoabsences
                                   PA.strategy = 'random') # Stratégie de sélection de pseudoabsences
  
  save(run.data, file = paste0("./", sp, "/run.data"))
  
  # Exploration du contenu
  run.data
  str(run.data, max.level = 3)
  
  

  # Illustration des pseudoabsences générées par biomod
  if(sp.list$pa.generation[i] == "biomod")
  {
    plot(baseline[["bio5"]])
    for (k in 1:ncol(run.data@PA))
    {
      tmp <- run.data@data.species
      tmp[is.na(tmp)] <- 2
      points(run.data@coord[run.data@PA[, k], ], pch = c(16, 1)[tmp[run.data@PA[, k]]], cex = .5, 
             col = brewer.pal(ncol(run.data@PA), "Paired")[k])
    }
  } 

 
  
  
  # 6. Calibration des modèles
  
  # k-fold cross validation : utiliser BIOMOD_cv
  # run.data.cv <- BIOMOD_cv(run.data,
                           # k = 3, repetition = 3)
  
  model.runs <- BIOMOD_Modeling(data = run.data, # Données initialisées par biomod
                                models =  c("GLM", # Liste des modèles à lancer
                                            "GAM",
                                            "ANN",
                                            "MARS",
                                            "GBM",
                                            "FDA",
                                            "RF"),
                                NbRunEval = 2, # Nombre de runs de cross-validation
                                DataSplit = 70, # % gardé pour la calibration dans la CV
                                # DataSplitTable = run.data.cv, # Pour la k-fold CV
                                VarImport = 3, # Nombre de runs de randomisation de variable importance
                                models.eval.meth = c("TSS", "ROC"), # Métriques d'évaluation (présence-absence)
                                do.full.models = FALSE, # Faire des modèles avec 100% des données ?
                                Prevalence = NULL, # Attribution de poids aux présences et absences pour forcer une prévalence
                                rescal.all.models = FALSE, # Faut-il rescaler les outputs des modèles dans le même range de valeurs ?
                                models.options = BIOMOD_ModelingOptions(
                                  MAXENT.Phillips = list(path_to_maxent.jar = "/home/pndb-cr/2-Galaxy/SDM/Cours_SDMs-master/data_cours")))

  # Sauvegarde des modèles calibrés sur le disque
  save(model.runs, file = paste0("./", sp, "/model.runs"))

  # Exploration du contenu
  str(model.runs)
  get_variables_importance(model.runs)
  model.runs@variables.importances@val
  varimp <- model.runs@variables.importances@val
  
  # 6. Préparation du modèle d'ensemble
  
  # Le modèle d'ensemble tel qu'effectué par biomod 
  # Il vaut mieux le construire soi-même
  # plutôt que d'utiliser biomod, par exemple pour évaluer les incertitudes
  # en profondeur ou pour se passer des limitations de cette fonction 
  # (i.e. ne permet pas d'inclure des répétitions hors biomod telles que des
  # pseudoabsences générées manuellement,
  # nombre de métriques limitées, ne calcule pas les écarts-types...)
  em.runs <- BIOMOD_EnsembleModeling(model.runs, # Modèles individuels calibrés
                                     chosen.models = 'all', # Filtrage des modèles pour l'EM
                                     em.by = 'all', # Comment construire l'EM ?
                                     # Voir la vignette EnsembleModelingAssembly dans l'aide
                                     eval.metric = 'TSS', # Quelle métrique utiliser pour filtrer les mauvais modèles 
                                     # ou pondérer la contribution des modèles dans l'EM ?
                                     eval.metric.quality.threshold = 0.6, # Quel seuil de filtration des mauvais modèles ? 
                                     models.eval.meth = c("TSS", "ROC"), # Quelles métriques utiliser pour évaluer l'EM ?
                                     prob.mean = TRUE, # Calculer la moyenne des probas
                                     prob.cv = TRUE, # Calculer le coefficient de variation
                                     prob.ci = TRUE, # Calculer l'intervalle de confiance autour de la moyenne
                                     prob.ci.alpha = 0.05, # Seuil de l'intervalle de confiance
                                     prob.median = TRUE, # Calculer la médiance
                                     committee.averaging = FALSE, # Faire du committee averaging (démocratie)
                                     prob.mean.weight = FALSE, # Calculer des moyennes pondérées
                                     prob.mean.weight.decay = 'proportional', # Méthode de pondération des modèles
                                     VarImport = 1) # Runs de variable importance pour l'EM
  
  # Sauvegarder l'EM sur le disque
  save(em.runs, file = paste0("./", sp, "/em.runs"))
  
  # # Exploration du contenu
  # str(em.runs, max.level = 5)
  # get_variables_importance(em.runs)
  # em.runs@em.models[[1]]@model_variables_importance
  # 
  
  
  for(j in listStacks)
  {
    # 7. Projection des modèles individuels
    
    cat(paste("---- ", Sys.time(), "Projection:", j, "----"))
    projection.stack <- readRDS(paste(initial.wd, "/data/", j, ".rds", sep = ""))
    projection.stack <- stack(projection.stack) # Au cas où le fichier est lu en tant que "brick" par raster
    # Dans ce cas, biomod n'aimera pas et ne fera pas la projection
    projection.runs <- BIOMOD_Projection(modeling.output = model.runs, # Modèles calibrés
                                         new.env = projection.stack[[cur.vars]], # Données environnementales sur lesquelles on projette les modèles
                                         proj.name = j, # Nom de la projection actuelle
                                         selected.models = 'all', # Modèles à projeter
                                         binary.meth = "TSS", # Avec quelle métrique faut-il transformer les probas en présence-absence ? 
                                         filtered.meth = NULL, # Avec quelle métrique appliquer un seuil en dessous duquel la proba est forcée à zéro ?
                                         build.clamping.mask = TRUE, # Le clamping mask illustre les zones où les prédictions sont en dehors des valeurs
                                         # utilisées lors de la calibration
                                         do.stack = T) # Stocker les projections dans un seul stack sur le disque ? Envisager de supprimer ce paramètre en cas
                                         # de fichiers trop gros
    save(projection.runs, file = paste("./", sp, "/proj_", j, "/projection.runs", sep = ""))
    
    # 8. Projection de l'EM
    proj.em <- BIOMOD_EnsembleForecasting(EM.output = em.runs, # Objet correspond à l'EM de biomod
                                          projection.output = projection.runs, # Objet issu des projections des modèles individuels
                                          binary.meth = "TSS", # Avec quelle métrique faut-il transformer les probas en présence-absence ? 
                                          do.stack = T) # Comme au-dessus
    save(proj.em, file = paste("./", sp, "/proj_", j, "/projection.em", sep = ""))
    cat(paste("---- ", Sys.time(), "Projection:", j, "finished ----\n"))
    print(warnings())
    
    # Vidage de la mémoire
    rm(list = c("projection.runs", "proj.em"))
    # Cette dernière ligne vide le répertoire temporaire de raster, 
    # d'expérience ça évite de se retrouver avec un disque plein sans s'en rendre compte...
    unlink(rasterOptions()$tmpdir, force = T)
  }
  unlink(dir(paste("./", sp, sep = ""), full.names = T)[grep("maxentWDtmp", dir(paste("./", sp, sep = "")))], force = T, recursive = T)
  setwd(initial.wd)
}

