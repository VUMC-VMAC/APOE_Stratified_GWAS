args <- commandArgs(TRUE)
filename <- args[1]  # Include path with file stem
output_location <- args[2]  # Output directory

library(qqman)
library(data.table)

options(bitmapType = 'cairo')

# Read in results
results <- fread(filename, header = TRUE, stringsAsFactors = FALSE)
names(results)[names(results) == "p-value"] <- "P"
names(results)[names(results) == "chromosome"] <- "CHR"
names(results)[names(results) == "position"] <- "BP"
names(results)[names(results) == "rs_number"] <- "SNP"

# Remove results with NA P values
if (sum(is.na(results$P)) > 0) {
  print(paste(sum(is.na(results$P)), "results with NA p-values. Removing..."))
  results <- results[!is.na(results$P), ]
}

# Remove results without CHR/BP info
if (sum(is.na(results$CHR)) > 0) {
  print(paste(sum(is.na(results$CHR)), "results without CHR/BP info. Removing to create the Manhattan plot..."))
  results <- results[!is.na(results$CHR), ]
}

# Get results for Manhattan plot
man_results <- results[, c("CHR", "BP", "SNP", "P")]

# Manhattan plot
output_filename <- file.path(output_location, paste0(basename(filename), ".manhattan.png"))
png(output_filename, width = 960, height = 480)
manhattan(man_results)
dev.off()

#calculate lambda
chisq <- qchisq(1-results$P,1)
lam <- median(chisq,na.rm=TRUE)/qchisq(0.5,1)

# QQ plot
output_filename_qq <- file.path(output_location, paste0(basename(filename), ".qq.png"))
png(output_filename_qq, width = 960, height = 480)
qq(results$P, sub = paste0("lambda = ", lam))
dev.off()
