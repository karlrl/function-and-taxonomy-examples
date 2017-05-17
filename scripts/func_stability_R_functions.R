# Load dependencies
library(vegan)


# Function to read in file, subset to particular samples, remove all rows with
# no nonzero values and get a vector of all sample pairwise spearman correlation
# coefficients and all Bray-Curtis Dissimilarity metrics
return_pairwise_coef <- function(filename, samples = NULL) {

  # Read in table
  table <- read.table(
    filename, header = T, sep = "\t", quote = "", stringsAsFactors = FALSE,
    row.names = 1)

  # Rename columns to remove ".fastq"
  colnames(table) <- sub(pattern = ".fastq", replacement = "", colnames(table))

  # Subset to specific samples if option set
  if(!is.null(samples)) {
    table <- table[, samples]
  }

  # Keep rows with with nonzero values in > 30% of samples:
  table_noZero <- table[rowSums(table > 0) > ceiling(ncol(table)*0.3), ]

  # Return list of spearman/brayCurtis metrics:
  return(pairwise_spearman_brayCurtis(table_noZero))
}


# Function to read in file, subset to particular samples, remove all rows with
# no nonzero values, then calculate and return [dis]similarity metrics
metaphlan2_cors <- function(filename, level, samples = NULL) {

  # Read in metaphlan2 output
  metaphlan2_out <- read.table(
    filename, header = T, sep = "\t" , stringsAsFactors = FALSE, quote = "")

  # Make vector that is all levels of interest concatenated
  metaphlan2_taxa <- c()
  for (j in 1:nrow(metaphlan2_out)) {
    metaphlan2_taxa <- c(
      metaphlan2_taxa, paste(metaphlan2_out[j, 1:level], collapse=";"))
  }

  # Remove first 8 columns (i.e. Kingdom:Strain)
  metaphlan2_rm <- metaphlan2_out[, -c(1:8)]

  # Set taxa column to be factor rather than character variable
  metaphlan2_rm$taxa <- as.factor(metaphlan2_taxa)

  # Sum each row that has the same value in taxa column
  metaphlan2_rm_summed <- aggregate(. ~ taxa, data = metaphlan2_rm, FUN = sum)

  # Set rownames equal to "taxa" column
  rownames(metaphlan2_rm_summed) <- metaphlan2_rm_summed$taxa

  # Remove "taxa" column
  metaphlan2_rm_summed <- metaphlan2_rm_summed[
    , -which(colnames(metaphlan2_rm_summed) == "taxa")]

  # Subset to specific samples if option set
  if(!is.null(samples)) {
    metaphlan2_rm_summed <- metaphlan2_rm_summed[, samples]
  }

  # Remove rows with with nonzero values in <= 50% of samples:
  metaphlan2_rm_summed_noZero <- metaphlan2_rm_summed[
    rowSums(metaphlan2_rm_summed > 0)
    > ceiling(ncol(metaphlan2_rm_summed)*0.3), ]

  # Return list of spearman/brayCurtis metrics:
  return(pairwise_spearman_brayCurtis(metaphlan2_rm_summed_noZero))
}


pairwise_spearman_brayCurtis <- function(input) {

  # Calculate Bray-Curtis dissimilarity
  input_bray = vegdist(t(input))

  # Calculate Jaccard dissimilarity
  input_jaccard = vegdist(t(input), method = "jaccard", binary = TRUE)

  # Calculate pairwise Spearman correlation coefficients between samples
  input_spearman <- cor(x = input, method = "spearman")

  # Set upper triagonal of matrix to be NA (including diagonal)
  input_spearman[upper.tri(input_spearman, diag = TRUE)] <- NA

  # Return list with two vectors: all pairwise spearman coefficients and all
  # Bray-Curtis Dissimilarity metrics
  return(list(
    "bray" = as.vector(input_bray),
    "jaccard" = as.vector(input_jaccard),
    "spearman" = as.vector(na.omit(as.vector(input_spearman)))
  ))
}
