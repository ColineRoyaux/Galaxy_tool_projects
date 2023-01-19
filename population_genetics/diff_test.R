#Rscript

################################################################################
##      Extract consensus sequence from aligned forward and reverse fasta     ##
################################################################################

#####Packages
library(adegenet, quietly = TRUE)
library(mmod, quietly = TRUE)

##Load arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("This tool needs at least one argument")
} else {

haploind <- args[1]
header <- as.logical(args[2])
haplo_col <- as.numeric(args[3]) #Column containing haplotypes
indivnames_col <- as.numeric(args[4]) #Column containing individuals
pop_col <- as.numeric(args[5]) #Column containing populations
}

##Steps
haploindiv_tab <- read.table(haploind, sep = "\t", header = header)
hap_tab <- as.data.frame(haploindiv_tab[, haplo_col])
row.names(hap_tab) <- haploindiv_tab[, indivnames_col]
colnames(hap_tab) <- "Locus"

#Convert dataset
hap_gen <- adegenet::df2genind(hap_tab, sep = "\t")

#Create population information
adegenet::pop(hap_gen) <- haploindiv_tab[, pop_col]

###Perform test
dif_t <- mmod::diff_test(hap_gen)

cat("Exact test of population genetic differentiation : \n\nPerformed on ", length(unique(haploindiv_tab[, indivnames_col])), "individuals of", length(unique(haploindiv_tab[, haplo_col])), "different haplotypes sampled in", length(unique(haploindiv_tab[, pop_col])), "populations.", "\n\np-value =", paste0(dif_t),"\n", file = "output.txt", append = TRUE)

if(dif_t <= 0.01){
    cat("\nSignificative test under 1%, haplotype composition is likely different between populations !\nHowever, check for potential bias in your data before making any conclusion.\n", file = "output.txt", append = TRUE)
}else if(dif_t <= 0.05){
    cat("\nSignificative test under 5%, haplotype composition is likely different between populations !\nHowever, check for potential bias in your data before making any conclusion.\n", file = "output.txt", append = TRUE)
}else{
    cat("\nTest isn't significative, no significative differences of haplotype composition between populations.\n", file = "output.txt", append = TRUE)
}

