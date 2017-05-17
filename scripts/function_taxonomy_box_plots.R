# Load in functions from Rscript
source(file = "scripts/func_stability_R_functions.R")

metadata <- read.csv(
  file = "metadata/HMASM_test_subset.csv",
  stringsAsFactors = FALSE
)

tongue_subset <- metadata[which(metadata$Body.Site == "tongue_dorsum"),
                          "SRS.ID"]
stool_subset <- metadata[which(metadata$Body.Site == "stool"),  "SRS.ID"]

# Get correleations for stool subset.
kos_stool_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_kos.spf", samples = stool_subset)
pathways_stool_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_pathways.spf", samples = stool_subset)
modules_stool_cor_vector <- return_pairwise_coef(
  filename ="tables/uniref100_to_modules.spf", samples = stool_subset)
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

# Get correlations for tongue subset
kos_tongue_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_kos.spf", samples = tongue_subset)
pathways_tongue_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_pathways.spf", samples = tongue_subset)
modules_tongue_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_modules.spf", samples = tongue_subset)
uniref100_tongue_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100.spf", samples = tongue_subset)
uniref90_tongue_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_uniref90.spf", samples = tongue_subset)
uniref50_tongue_cor_vector <- return_pairwise_coef(
  filename = "tables/uniref100_to_uniref50.spf", samples = tongue_subset)

strains_tongue_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 8, samples = tongue_subset)
species_tongue_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 7, samples = tongue_subset)
genus_tongue_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 6, samples = tongue_subset)
family_tongue_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 5, samples = tongue_subset)
order_tongue_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 4, samples = tongue_subset)
class_tongue_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 3, samples = tongue_subset)
phylum_tongue_cor_vector <- metaphlan2_cors(
  filename = "tables/metaphlan-taxonomy.spf", level = 2, samples = tongue_subset)

spearman_combined_stool_df <- data.frame(
  "UniRef100" = uniref100_stool_cor_vector$spearman,
  "UniRef90" = uniref90_stool_cor_vector$spearman,
  "UniRef50" = uniref50_stool_cor_vector$spearman,
  "KOs" = kos_stool_cor_vector$spearman,
  "Modules" = modules_stool_cor_vector$spearman,
  "Pathways" = pathways_stool_cor_vector$spearman,
  "Strains" = strains_stool_cor_vector$spearman,
  "Species" = species_stool_cor_vector$spearman,
  "Genus" = genus_stool_cor_vector$spearman,
  "Family" = family_stool_cor_vector$spearman,
  "Order" = order_stool_cor_vector$spearman,
  "Class" = class_stool_cor_vector$spearman,
  "Phylum" = phylum_stool_cor_vector$spearman
)

brayCurtis_combined_stool_df <- data.frame(
  "UniRef100" = uniref100_stool_cor_vector$bray,
  "UniRef90" = uniref90_stool_cor_vector$bray,
  "UniRef50" = uniref50_stool_cor_vector$bray,
  "KOs" = kos_stool_cor_vector$bray,
  "Modules" = modules_stool_cor_vector$bray,
  "Pathways" = pathways_stool_cor_vector$bray,
  "Strains" = strains_stool_cor_vector$bray,
  "Species" = species_stool_cor_vector$bray,
  "Genus" = genus_stool_cor_vector$bray,
  "Family" = family_stool_cor_vector$bray,
  "Order" = order_stool_cor_vector$bray,
  "Class" = class_stool_cor_vector$bray,
  "Phylum" = phylum_stool_cor_vector$bray
)

jaccard_combined_stool_df <- data.frame(
  "UniRef100" = uniref100_stool_cor_vector$jaccard,
  "UniRef90" = uniref90_stool_cor_vector$jaccard,
  "UniRef50" = uniref50_stool_cor_vector$jaccard,
  "KOs" = kos_stool_cor_vector$jaccard,
  "Modules" = modules_stool_cor_vector$jaccard,
  "Pathways" = pathways_stool_cor_vector$jaccard,
  "Strains" = strains_stool_cor_vector$jaccard,
  "Species" = species_stool_cor_vector$jaccard,
  "Genus" = genus_stool_cor_vector$jaccard,
  "Family" = family_stool_cor_vector$jaccard,
  "Order" = order_stool_cor_vector$jaccard,
  "Class" = class_stool_cor_vector$jaccard,
  "Phylum" = phylum_stool_cor_vector$jaccard
)

spearman_combined_tongue_df <- data.frame(
  "UniRef100" = uniref100_tongue_cor_vector$spearman,
  "UniRef90" = uniref90_tongue_cor_vector$spearman,
  "UniRef50" = uniref50_tongue_cor_vector$spearman,
  "KOs" = kos_tongue_cor_vector$spearman,
  "Modules" = modules_tongue_cor_vector$spearman,
  "Pathways" = pathways_tongue_cor_vector$spearman,
  "Strains" = strains_tongue_cor_vector$spearman,
  "Species" = species_tongue_cor_vector$spearman,
  "Genus" = genus_tongue_cor_vector$spearman,
  "Family" = family_tongue_cor_vector$spearman,
  "Order" = order_tongue_cor_vector$spearman,
  "Class" = class_tongue_cor_vector$spearman,
  "Phylum" = phylum_tongue_cor_vector$spearman
)

brayCurtis_combined_tongue_df <- data.frame(
  "UniRef100" = uniref100_tongue_cor_vector$bray,
  "UniRef90" = uniref90_tongue_cor_vector$bray,
  "UniRef50" = uniref50_tongue_cor_vector$bray,
  "KOs" = kos_tongue_cor_vector$bray,
  "Modules" = modules_tongue_cor_vector$bray,
  "Pathways" = pathways_tongue_cor_vector$bray,
  "Strains" = strains_tongue_cor_vector$bray,
  "Species" = species_tongue_cor_vector$bray,
  "Genus" = genus_tongue_cor_vector$bray,
  "Family" = family_tongue_cor_vector$bray,
  "Order" = order_tongue_cor_vector$bray,
  "Class" = class_tongue_cor_vector$bray,
  "Phylum" = phylum_tongue_cor_vector$bray
)

jaccard_combined_tongue_df <- data.frame(
  "UniRef100" = uniref100_tongue_cor_vector$jaccard,
  "UniRef90" = uniref90_tongue_cor_vector$jaccard,
  "UniRef50" = uniref50_tongue_cor_vector$jaccard,
  "KOs" = kos_tongue_cor_vector$jaccard,
  "Modules" = modules_tongue_cor_vector$jaccard,
  "Pathways" = pathways_tongue_cor_vector$jaccard,
  "Strains" = strains_tongue_cor_vector$jaccard,
  "Species" = species_tongue_cor_vector$jaccard,
  "Genus" = genus_tongue_cor_vector$jaccard,
  "Family" = family_tongue_cor_vector$jaccard,
  "Order" = order_tongue_cor_vector$jaccard,
  "Class" = class_tongue_cor_vector$jaccard,
  "Phylum" = phylum_tongue_cor_vector$jaccard
)

label_col <- c(
  replicate(6, "grey"),  # Functional "levels"
  replicate(8, "white")  # Taxonomic levels
)

boxplot(
  spearman_combined_stool_df,
  ylim = c(-0.1, 1),
  ylab = "Spearman's rho",
  outline = TRUE,
  main = "Stool - Pairwise Spearman",
  col = label_col,
  las = 2
)
boxplot(
  brayCurtis_combined_stool_df,
  ylim = c(0, 1),
  ylab = "Bray-Curtis Dissimilarity",
  outline = TRUE,
  main = "Stool - Pairwise Bray-Curtis",
  col = label_col,
  las = 2
)
boxplot(
  jaccard_combined_stool_df,
  ylim = c(0, 1),
  ylab = "Jaccard Index",
  outline = TRUE,
  main = "Stool - Jaccard",
  col = label_col,
  las = 2
)
boxplot(
  spearman_combined_tongue_df,
  ylim = c(-0.1, 1),
  ylab = "Spearman's rho",
  outline = TRUE,
  main = "Tongue - Pairwise Spearman",
  col = label_col, las = 2
)
boxplot(
  brayCurtis_combined_tongue_df,
  ylim = c(0, 1),
  ylab = "Bray-Curtis Dissimilarity",
  outline = TRUE,
  main = "Tongue - Pairwise Bray-Curtis",
  col = label_col,
  las = 2
)
boxplot(
  jaccard_combined_tongue_df,
  ylim = c(0, 1),
  ylab = "Jaccard Index",
  outline = TRUE,
  main = "Tongue - Jaccard",
  col = label_col,
  las = 2
)

