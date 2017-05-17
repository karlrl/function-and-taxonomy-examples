calc_variance_from_aggregate <- function(mean_val, freq, val) {
  
  diff_sq_sum <- 0
  for(i in val){
    diff_sq_sum <- diff_sq_sum + freq[i]*(i-mean_val)**2
  }
  
  return(diff_sq_sum/(sum(freq)-1))
}

num_taxa_barplot2_input <- function(input_file) {

  input_table <- data.frame(t(read.table(input_file, header=T ,sep="\t", 
                                       row.names=1)))
                                        
  # Remove "0" counts:
  input_table <- input_table[-1,]
  
  # Remove last (empty) line:
  input_table <- input_table[-nrow(input_table),]
  
  # Get bins in numeric vector
  bin_val <- as.numeric(sub("^X", "", rownames(input_table)))
  
  # Intitialize empty vectors to be returned.
  func_means <- c()
  func_lower_ci <- c()
  func_upper_ci <- c()
  func_var <- c()
  func_sem <- c()
  func_names <- c()
  
  for (func_category in colnames(input_table)) {
    
    cat_mean <- sum(input_table[,func_category]*bin_val)/sum(input_table[,func_category])
    
    cat_variance <- calc_variance_from_aggregate(mean_val=cat_mean, 
                                                 freq=input_table[,func_category],
                                                 val=bin_val)
    
    cat_sem <- sqrt(cat_variance)/sqrt(sum(input_table[,func_category]))
    
    func_means <- c(func_means, cat_mean)
    func_lower_ci <- c(func_lower_ci, cat_mean-2*cat_sem)
    func_upper_ci <- c(func_upper_ci, cat_mean+2*cat_sem)
    func_var <- c(func_var, cat_variance)
    func_sem <- c(func_sem, cat_sem)
    func_names <- c(func_names, func_category)
  }
 
  func_barplot2_input <- list(means=func_means, lower_ci = func_lower_ci,
                       upper_ci = func_upper_ci, var=func_var, sem=func_sem,
                       names = func_names)
  
  return(func_barplot2_input)

}
  