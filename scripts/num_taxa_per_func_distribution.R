library(gplots)

source("scripts/plot_mean_num_taxa_R_functions.R")

species_input <- num_taxa_barplot2_input("tables/species_function_counts.txt")

superkingdom_input <- num_taxa_barplot2_input("tables/superkingdom_function_counts.txt")

cat_names <- c("UniRef100", "UniRef90", "UniRef50", "KEGG\nOrthologs",
               "KEGG\nModules", "KEGG\nPathways")

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

# kingdom_input <- num_taxa_barplot2_input("kingdom_function_counts.txt")
# barplot2(height = kingdom_input$means, names.arg = cat_names, plot.ci = TRUE,
#          ci.l = kingdom_input$lower_ci, ci.u = kingdom_input$upper_ci,
#          ylim=c(0,2), las=3,
#          ylab="Mean number of kingdoms that possess function")


### Converting back to distribution to double-check calculations:
# tmp <- data.frame(t(read.table("species_function_counts.txt", header=T ,sep="\t",
#                                row.names=1)))
# tmp <- tmp[-1,]
# tmp <- tmp[-nrow(tmp),]
# xlab_val <- as.numeric(sub("^X", "", rownames(tmp)))
# ko_dist <- c()
# for(i in xlab_val){
#   if(as.numeric(tmp[i, "ko"]) > 0) {
#     ko_dist <- c(ko_dist, rep.int(i,tmp[i, "ko"]))
#   }
# }

dev.off()
