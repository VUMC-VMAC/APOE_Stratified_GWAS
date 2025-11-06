args <- commandArgs(TRUE)
filename <- args[1]  # Include path with file stem
snp_info <- args[2]  # BIM file if these are meta-analysis results

# Set plot directory
plot_dir <- "/data/h_vmac/contra2/01_APOE4_Strat/05_MAGMA/plots/"

# Load necessary libraries
library(qqman)
library(data.table)

# Read in results
results <- fread(filename, header=TRUE, stringsAsFactors=FALSE)
names(results)[names(results) == "p-value"] <- "P"
names(results)[names(results) == "chromosome"] <- "CHR"
names(results)[names(results) == "START"] <- "BP"
names(results)[names(results) == "rs_number"] <- "SNP"

# Check if chromosome and position columns are present; if not, add from SNP info
if(!("CHR" %in% names(results)) || !("BP" %in% names(results))){
  print("Chromosome and position not present in results dataframe. Pulling in now...")
  snps <- fread(snp_info, header=FALSE, stringsAsFactors=FALSE)
  snps <- snps[,c(1,2,4)]
  names(snps) <- c("CHR", "SNP", "BP")
  snps <- snps[snps$BP %in% results$BP,]
  results <- merge(results, snps, by="BP", all.x=TRUE)
}

# Remove results with NA P values
if(sum(is.na(results$P)) > 0){
  print(paste(sum(is.na(results$P)), "results with NA p values. Removing..."))
  results <- results[!is.na(results$P),]
}

# Calculate lambda
chisq <- qchisq(1 - results$P, 1)
lam <- median(chisq, na.rm=TRUE) / qchisq(0.5, 1)

# Generate Q-Q plot and save to plot directory
print("Making Q-Q plot")
qq_plot_path <- paste0(plot_dir, basename(filename), ".qq.png")
png(qq_plot_path, type="cairo")
qq(results$P, sub=paste0("lambda=", lam))
dev.off()

# Remove results without CHR/BP info for Manhattan plot
if(sum(is.na(results$CHR)) > 0){
  print(paste(sum(is.na(results$CHR)), "results without CHR/BP info. Removing to create the Manhattan plot..."))
  results <- results[!is.na(results$CHR),]
}

# Prepare data for Manhattan plot
man_results <- results[,c("CHR", "BP", "P")]

# Generate Manhattan plot and save to plot directory
print("Making Manhattan plot")
manhattan_plot_path <- paste0(plot_dir, basename(filename), ".manhattan.png")
png(manhattan_plot_path, width=480, height=480, type="cairo")
manhattan(man_results, ylim=c(0,10), cex.axis=1.2, genomewideline=FALSE)
dev.off()
