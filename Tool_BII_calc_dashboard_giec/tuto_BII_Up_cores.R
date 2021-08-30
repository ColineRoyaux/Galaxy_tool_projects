suppressMessages(library(dplyr)) # for easy data manipulation
suppressMessages(library(tidyr)) # ditto
suppressMessages(library(magrittr)) # for piping
suppressMessages(library(lme4)) # for mixed effects models
suppressMessages(library(car)) # for logit transformation with adjustment
suppressMessages(library(raster)) # for working with raster data
suppressMessages(library(geosphere)) # calculating geographic distance between sites
suppressMessages(library(foreach)) # running loops
suppressMessages(library(doParallel)) # running loops in parallel
suppressMessages(library(ggplot2)) # nice plots

args = commandArgs(trailingOnly=TRUE)

## Get land use data
unzip(args[2], exdir = "land_use")
#stop(list.dirs(getwd(), recursive = TRUE))

files_zip <- grep(".zip$", list.files("land_use", recursive = TRUE, full.names = TRUE), value = TRUE)

if(!is.null(files_zip)) {
    l <- lapply(files_zip, 
                function(x) {
                    unzip(x, exdir = "land_use/layers")
                })
}

files_G <- grep(".bil$", list.files("land_use/layers", recursive = TRUE, full.names = TRUE), value = TRUE)

rast_from_files <- lapply(files_G, #Try to create a raster from each listed files
                          function(x) {
                            r <- tryCatch(raster(x), error = function(e) {}) #Need to clear the  stdout/stderr
                            return(r)
                          })

rast_names <- lapply(rast_from_files,
                     function(x) {
                       return(names(x))
                     })

names(rast_from_files) <- substr(rast_names, 1, 3)

dir.create("land_use/plots")

rast_red_f <- lapply(rast_from_files,
                     function (x) {
                       ## Global rasters are too big we need to get a smaller extent to have a data frame < 500Mb: Europe extent(-11, 50, 35, 60), paleartic west extent(-15, 50, 40, 70), France extent(-5, 10, 41, 52), retained extent(-11, 42, 36, 60), Americas extent(-160, -20, -60, 65)
                       rast_red <- crop(x, extent(-11, 50, 35, 70))
                       rast_df<- as.data.frame(rast_red, xy=TRUE)
                       names(rast_df) <- c("x", "y", "layer")
                       p <- ggplot(data = rast_df) +
                            geom_raster(mapping=aes(x=x, y=y, fill=layer)) +
                            scale_fill_gradientn(colours= rev(terrain.colors(10)), name = paste0(substr(names(x), 1, 3), " landuse proportion"), na.value = "sky blue") +
                            theme(panel.background = element_blank(), panel.grid = element_line(colour = TRUE))
                      ggsave(paste0("land_use/plots/plot_init_", substr(names(x), 1, 3),".png"), p, width = 25, height = 15, units = "cm")
                      return(rast_red)
       }
)


# read in the data
diversity <- readRDS(args[1]) %>%
  
  # now let's filter out just the data for the Americas
  filter(UN_region == "Europe")

glimpse(diversity)


table(diversity$Predominant_land_use, diversity$Use_intensity)

diversity <- diversity %>%
  # make a level of Primary minimal. Everything else gets the coarse land use
  mutate(
    LandUse = ifelse(Predominant_land_use == "Primary vegetation" & Use_intensity == "Minimal use",
                     "Primary minimal",
                     paste(Predominant_land_use)),
    
    # collapse the secondary vegetation classes together
    LandUse = ifelse(grepl("secondary", tolower(LandUse)),
                     "Secondary vegetation",
                     paste(LandUse)),
    
    # change cannot decide into NA
    LandUse = ifelse(Predominant_land_use == "Cannot decide",
                     NA, 
                     paste(LandUse)),
    
    # relevel the factor so that Primary minimal is the first level (so that it is the intercept term in models)
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary minimal")
  )


abundance_data <- diversity %>%
  
  # pull out just the abundance measures
  filter(Diversity_metric_type == "Abundance") %>%
  
  # group by SSBS (each unique value corresponds to a unique site)
  group_by(SSBS) %>%
  
  # now add up all the abundance measurements within each site
  mutate(TotalAbundance = sum(Effort_corrected_measurement)) %>%
  
  # ungroup
  ungroup() %>%
  
  # pull out unique sites
  distinct(SSBS, .keep_all = TRUE) %>%
  
  # now group by Study ID
  group_by(SS) %>%
  
  # pull out the maximum abundance for each study
  mutate(MaxAbundance = max(TotalAbundance)) %>%
  
  # ungroup
  ungroup() %>%
  
  # now rescale total abundance, so that within each study, abundance varies from 0 to 1.
  mutate(RescaledAbundance = TotalAbundance/MaxAbundance)

cd_data_input <- diversity %>%
  
  # drop any rows with unknown LandUse
  filter(!is.na(LandUse)) %>%
  
  # pull out only the abundance data
  filter(Diversity_metric_type == "Abundance") %>%
  
  # group by Study
  group_by(SS) %>%
  
  # calculate the number of unique sampling efforts within that study
  mutate(n_sample_effort = n_distinct(Sampling_effort)) %>%
  
  # calculate the number of unique species sampled in that study
  mutate(n_species = n_distinct(Taxon_name_entered)) %>%
  
  # check if there are any Primary minimal sites in the dataset
  mutate(n_primin_records = sum(LandUse == "Primary minimal")) %>%
  
  # ungroup
  ungroup() %>%
  
  # now keep only the studies with one unique sampling effort
  filter(n_sample_effort == 1) %>%
  
  # and keep only studies with more than one species 
  # as these studies clearly aren't looking at assemblage-level diversity
  filter(n_species > 1) %>%
  
  # and keep only studies with at least some Primary minimal data
  filter(n_primin_records > 0) %>%
  
  # drop empty factor levels
  droplevels()

getJacAbSym <- function(s1, s2, data){

  # get the list of species that are present in site 1 (i.e., their abundance was greater than 0)
  s1species <- data %>%
    
    # filter out the SSBS that matches s1
    filter(SSBS == s1) %>%
    
    # filter out the species where the Measurement (abundance) is greater than 0
    filter(Measurement > 0) %>%
    
    # get the unique species from this dataset
    distinct(Taxon_name_entered) %>%
    
    # pull out the column into a vector
    pull
  
  # for site 2, get the total abundance of species that are also present in site 1
  
  s2abundance_s1species <- data %>%
    
    # filter out the SSBS that matches s2
    filter(SSBS == s2) %>%
    
    # filter out the species that are also present in site 1
    filter(Taxon_name_entered %in% s1species) %>%
    
    # pull out the Measurement into a vector
    pull(Measurement) %>%
    
    # calculate the sum
    sum()
  
  # calculate the total abundance of all species in site 2
  s2_sum <- data %>%
    
    # filter out the SSBS that matches s2
    filter(SSBS == s2) %>%
    
    # pull out the measurement column (the abundance)
    pull(Measurement) %>%
    
    # calculate the sum
    sum() 
  
  
  # Now calculate the compositional similarity
  # this is the number of individuals of species also found in s1, divided by the total abundance in s2 
  # so that equates to the proportion of individuals in s2 that are of species also found in s1
  
  sor <- s2abundance_s1species / s2_sum

  
  # if there are no taxa in common, then sor = 0
  # if abundances of all taxa are zero, then similarity becomes NaN.
  return(sor)
  
  }


# get a vector of each study to loop over
studies <- distinct(cd_data_input, SS) %>%
  pull()

cl <- detectCores() %>% -1 %>% makeCluster

registerDoParallel(cores = cl)

# If you're not familiar with loops (or with foreach loops):
# I'm going to loop over every element (s) in studies and combine the results of each loop by rbinding them into one large dataset. Since we're using functions from different packages within this loop, we need to specify them (if you don't do the loop in parallel, this isn't necessary)
cd_data <- foreach(s = studies, 
                   .combine = rbind,
                   .packages = c("dplyr", "magrittr", "geosphere")) %dopar% {
  
  # filter out the given study
  data_ss <- filter(cd_data_input, SS == s)
  
  # pull out the SSBS and LandUse information (we'll use this later to assign a land use contrast to each pair of site
  site_data <- data_ss %>%
    dplyr::select(SSBS, LandUse) %>%
    distinct(SSBS, .keep_all = TRUE)
  
  # pull out the sites that are Primary minimal (we only want to use comparisons with the baseline)
  baseline_sites <- site_data %>%
    filter(LandUse == "Primary minimal") %>%
    pull(SSBS)
  
  # pull out all the sites
  site_list <- site_data %>%
    pull(SSBS)
  

  # get all site x site comparisons for this study
  site_comparisons <- expand.grid(baseline_sites, site_list) %>%
    
    # rename the columns so they will be what the compositional similarity function expects for ease
    rename(s1 = Var1, s2 = Var2) %>%
    
    # remove the comparisons where the same site is being compared to itself
    filter(s1 != s2)
    
  
  # apply the compositional similarity function over each site combination in the dataset
  sor <- apply(site_comparisons, 1, function(y) getJacAbSym(data = data_ss, s1 = y['s1'], s2 = y['s2']))
  
  # calculate the geographic distance between sites
  # first pull out the lat and longs for each site combination
  s1LatLong <- as.matrix(data_ss[match(site_comparisons$s1, data_ss$SSBS), c('Longitude','Latitude')])
  s2LatLong <- as.matrix(data_ss[match(site_comparisons$s2, data_ss$SSBS), c('Longitude','Latitude')])

  # then calculate the distance between sites
  dist <- distHaversine(s1LatLong, s2LatLong)

  # pull out the land-use contrast for those site combinations
  Contrast <- paste(site_data$LandUse[match(site_comparisons$s1, site_data$SSBS)],
                    site_data$LandUse[match(site_comparisons$s2, site_data$SSBS)], 
                    sep = "-")
  
  # put all the information into a single dataframe
  
  study_results <- data.frame(site_comparisons,
                              sor,
                              dist,
                              Contrast,
                              SS = s,
                              stringsAsFactors = TRUE)

  
  
}

# stop running things in parallel
registerDoSEQ()


# run a simple model
ab_m <- lmer(sqrt(RescaledAbundance) ~ LandUse + (1|SS) + (1|SSB), data = abundance_data)
capture.output(summary(ab_m), file = "abundance_model.txt")


# there is some data manipulation we want to do before modelling
cd_data <- cd_data %>%
  
  # Firstly, we only care about comparisons where Primary minimal is the first site
  # so pull the contrast apart
  separate(Contrast, c("s1_LandUse", "s2_LandUse"), sep = "-", remove = FALSE) %>%
  
  # filter sites where s1_LandUse is the baseline site
  filter(s1_LandUse == "Primary minimal") %>%
  
  # logit transform the compositional similarity
  mutate(logitCS = logit(sor, adjust = 0.001, percents = FALSE)) %>%
  
  # log10 transform the geographic distance between sites
  mutate(log10geo = log10(dist + 1)) %>%
  
  # make primary minimal-primary minimal the baseline again
  mutate(Contrast = factor(Contrast), 
         Contrast = relevel(Contrast, ref = "Primary minimal-Primary minimal"))


# Model compositional similarity as a function of the land-use contrast and the geographic distance between sites
cd_m <- lmer(logitCS ~ Contrast + log10geo + (1|SS), data = cd_data)
capture.output(summary(cd_m), file = "similarity_model.txt")


# let's start with the abundance model

# set up a dataframe with all the levels you want to predict diversity for
# so all the land-use classes in your model must be in here
newdata_ab <- data.frame(LandUse = levels(abundance_data$LandUse)) %>%
  
  # now calculate the predicted diversity for each of these land-use levels
  # setting re.form = NA means random effect variance is ignored
  # then square the predictions (because we modelled the square root of abundance, so we have to back-transform it to get the real predicted values)
  mutate(ab_m_preds = predict(ab_m, ., re.form = NA) ^ 2)

#stop(newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Primary minimal'])

# now the compositional similarity model

# set up a function to calculate the inverse logit of the adjusted logit function we used
# where f is the value to be back-transformed and a is the adjustment value used for the transformation
inv_logit <- function(f, a){
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

# once again, set up the dataframe with all the levels you want to predict diversity for
# because we had an extra fixed effect in this model (log10geo), we also have to set a baseline level for this
# we'll set it to 0 because we're interested in the compositional similarity when distance between sites is 0 (i.e., when distance-based turnover is discounted, what is the turnover because of land use?)
newdata_cd <- data.frame(Contrast = levels(cd_data$Contrast),
                      log10geo = 0) %>%
  mutate(cd_m_preds = predict(cd_m, ., re.form = NA) %>%
           inv_logit(a = 0.001))



# generate a dataframe with random numbers for the cell values for each land-use class
#lus <- data.frame(pri_min = rnorm(25, mean = 50, sd = 25),
 #                 pri = rnorm(25, mean = 100, sd = 25),
  #                plant = rnorm(25, mean = 100, sd = 25),
   #               sec = rnorm(25, mean = 300, sd = 25),
    #              crop = rnorm(25, mean = 1000, sd = 25),
     #             pas = rnorm(25, mean = 400, sd = 25),
      #            urb = rnorm(25, mean = 50, sd = 25)
       #           )

# let's artificially make the first cell dominated by urban land and the last cell dominated by minimally-used primary vegetation
#lus$urb[1] <- 2000
#lus$pri_min[25] <- 2000

#lus <- lus %>%
  # calculate the row totals
 # mutate(tot = rowSums(.)) %>%
  
  # now, for each land use, divide the value by the rowsum
  # this will give us the proportion of each land use in each cell
  #transmute_at(1:7, list(~ ./tot))

# double check that the proportions of each land use sum to 1 (accounting for rounding errors)
#cat("Does the proportions of each land use sum up to 1 ?", all(zapsmall(rowSums(lus)) == 1), "If FALSE, there must be a problem")


# for each column of lus (i.e., each land use)
#for(i in 1:ncol(lus)){
  
  # take the column and turn it into a 5 x 5 matrix
 # ras <- matrix(lus[ , i], nrow = 5, ncol = 5) %>%
    # turn that into a raster
  #  raster()
  
  # come up with a name for the object to hold that raster
  #nm <- paste(names(lus)[i], "raster", sep = "_")
  
  # and assign the raser to that name
  #assign(x = nm, value = ras)
  
#}

#stop(pri_min_raster)
#plot(urb_raster)
#png(paste0(getwd(), "raster_urban_land.png"), width = 550, height = 550, res = 300)

ab_raster <- (newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Primary vegetation'] * rast_red_f$PRI + 
  newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Secondary vegetation'] * rast_red_f$SEC +
  newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Cropland'] * rast_red_f$CRP +
  newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Pasture'] * rast_red_f$PAS +
  newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Urban'] * rast_red_f$URB) /
  
  # divide by the reference value
  newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Primary minimal']


#ab_raster <- (newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Primary minimal'] * pri_min_raster + 
 # newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Primary vegetation'] * pri_raster + 
  #newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Plantation forest'] * plant_raster +
#  newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Secondary vegetation'] * sec_raster +
 # newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Cropland'] * crop_raster +
  #newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Pasture'] * pas_raster +
  #newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Urban'] * urb_raster) /
  
  # divide by the reference value
  #newdata_ab$ab_m_preds[newdata_ab$LandUse == 'Primary minimal']

cd_raster <- (newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Primary vegetation'] * rast_red_f$PRI + 
  newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Secondary vegetation'] * rast_red_f$SEC +
  newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Cropland'] * rast_red_f$CRP +
  newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Pasture'] * rast_red_f$PAS +
  newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Urban'] * rast_red_f$URB) /
  
  # divide by the reference value
  newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Primary minimal']

#cd_raster <- (newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Primary minimal'] * pri_min_raster + 
 # newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Primary vegetation'] * pri_raster + 
  #newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Plantation forest'] * plant_raster + 
#  newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Secondary vegetation'] * sec_raster +
 # newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Cropland'] * crop_raster +
  #newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Pasture'] * pas_raster +
#  newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Urban'] * urb_raster) /
  
  # divide by the reference value
 # newdata_cd$cd_m_preds[newdata_cd$Contrast == 'Primary minimal-Primary minimal']

bii <- ab_raster * cd_raster
#plot(bii * 100)
#png(paste0(getwd(), "BII.png"), width = 550, height = 550, res = 300)

biib <- bii * 100
#plot(biib)

bii_crop <- crop(biib, extent(-11, 50, 35, 70))

bii_rast <- as.data.frame(bii_crop, xy=TRUE)

pl <- ggplot(data = bii_rast)+
      geom_raster(mapping=aes(x=x, y=y, fill=layer))+
      scale_fill_gradientn(colours= rev(terrain.colors(10)), name='BII', na.value = "sky blue") +
      theme(panel.background = element_blank(), panel.grid = element_line(colour = TRUE))
ggsave("BII.png", pl, width = 25, height = 15, units = "cm")
#utils::sessionInfo()


