library(gplots)

# Load in accessory functions from other Rscript.
source("scripts/plot_mean_num_taxa_R_functions.R")

# Load in input files and get mean and standard error for each functional type.
species_input <- num_taxa_barplot2_input("tables/species_function_counts.txt")
superkingdom_input <- num_taxa_barplot2_input("tables/superkingdom_function_counts.txt")

# Specify order of func types.
cat_names <- c("UniRef100", "UniRef90", "UniRef50", "KEGG\nOrthologs",
               "KEGG\nModules", "KEGG\nPathways")

# Plot barplots and write to pdf as two panels.
pdf("figures/num_taxa_per_func_distribution.pdf", width = 11, height = 8.5)
par(mfrow=c(1,2))

barplot2(height = species_input$means, names.arg = cat_names, plot.ci = TRUE,
         ci.l = species_input$lower_ci, ci.u = species_input$upper_ci,
         ylim=c(0,2000), las=3,
         ylab="Mean number of species that possess function")

barplot2(height = superkingdom_input$means, names.arg = cat_names, plot.ci = TRUE,
         ci.l = superkingdom_input$lower_ci, ci.u = superkingdom_input$upper_ci,
         ylim=c(0,3), las=3,
         ylab="Mean number of superkingdoms that possess function")

dev.off()
