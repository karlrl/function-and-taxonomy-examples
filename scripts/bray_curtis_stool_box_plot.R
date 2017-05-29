# Load in functions from Rscript
source(file = "scripts/func_stability_R_functions.R")

metadata <- read.csv(
  file = "metadata/HMASM_test_subset.csv",
  stringsAsFactors = FALSE
)

stool_subset <- metadata[which(metadata$Body.Site == "stool"),  "SRS.ID"]

# Get correleations for stool subset.
kos_stool_cor_vector <- return_pairwise_coef(
  filename = "tables/humann_kos.spf", samples = stool_subset)
pathways_stool_cor_vector <- return_pairwise_coef(
  filename = "tables/humann_pathways.spf", samples = stool_subset)
modules_stool_cor_vector <- return_pairwise_coef(
  filename ="tables/humann_modules.spf", samples = stool_subset)
uniref100_stool_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100.spf", samples = stool_subset)
uniref90_stool_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_uniref90.spf", samples = stool_subset)
uniref50_stool_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_uniref50.spf", samples = stool_subset)

strains_stool_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 8, samples = stool_subset)
species_stool_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 7, samples = stool_subset)
genus_stool_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 6, samples = stool_subset)
family_stool_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 5, samples = stool_subset)
order_stool_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 4, samples = stool_subset)
class_stool_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 3, samples = stool_subset)
phylum_stool_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 2, samples = stool_subset)

brayCurtis_combined_stool_df <- data.frame(
  "UniRef100" = uniref100_stool_cor_vector$bray,
  "UniRef90" = uniref90_stool_cor_vector$bray,
  "UniRef50" = uniref50_stool_cor_vector$bray,
  "KEGG Orthologs" = kos_stool_cor_vector$bray,
  "KEGG Modules" = modules_stool_cor_vector$bray,
  "KEGG Pathways" = pathways_stool_cor_vector$bray,
  "Strains" = strains_stool_cor_vector$bray,
  "Species" = species_stool_cor_vector$bray,
  "Genus" = genus_stool_cor_vector$bray,
  "Family" = family_stool_cor_vector$bray,
  "Order" = order_stool_cor_vector$bray,
  "Class" = class_stool_cor_vector$bray,
  "Phylum" = phylum_stool_cor_vector$bray
)

label_col <- c(
  replicate(6, "grey"),  # Functional "levels"
  replicate(8, "white")  # Taxonomic levels
)

labels <- c(
  "UniRef100",
  "UniRef90",
  "UniRef50",
  "KEGG\nOrthologs",
  "KEGG\nModules",
  "KEGG\nPathways",
  "Strains",
  "Species",
  "Genus",
  "Family",
  "Order",
  "Class",
  "Phylum"
)

# 5.5 x 8 is good for pasting in Word (with cex=1)
par(cex=1.5)
pdf("bray_curtis_stool.pdf", height = 11, width = 8.5)

boxplot(
  brayCurtis_combined_stool_df,
  ylim = c(0, 1),
  ylab = "Bray-Curtis Dissimilarity",
  names = labels,
  outline = TRUE,
  col = label_col,
  las = 2
)
legend(
  "topright",
  inset = 0.02,
  c("Function", "Taxonomy"),
  bty = "o",
  fill=c("grey", "white")
)
