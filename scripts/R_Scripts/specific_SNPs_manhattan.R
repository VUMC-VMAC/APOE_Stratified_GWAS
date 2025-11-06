args <- commandArgs(TRUE)
filename <- "/data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/results/MetaCrossAnc_CrossAnc_EXF_APOE4pos.in_2of2Ancestry.out"
snp_info <- "/data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files/AllRaces/AllRaces_combined.bim"

 #/data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files/NHW/NHW_combined.bim"

# Bim: /data/h_vmac/contra2/01_APOE4_Strat/00_Genetic_Files/AllRaces/AllRaces_combined.bim
# /data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/results/MetaCrossAnc_CrossAnc_EXF_APOE4neg.in_2of2Ancestry.out
# /data/h_vmac/contra2/01_APOE4_Strat/04_CrossAnc/results/MetaCrossAnc_CrossAnc_EXF_APOE4pos.in_2of2Ancestry.out

#Run this
# /data/h_vmac/contra2/01_APOE4_Strat/scripts/R_Scripts/specific_SNPs_manhattan.R

library(qqman)
library(data.table)

# Define output directory
output_dir <- "/data/h_vmac/contra2/01_APOE4_Strat/scripts/R_Scripts/"

# Read in results
results <- fread(filename, header=TRUE, stringsAsFactors = FALSE)
names(results)[names(results) == "p-value"] <- "P"
names(results)[names(results) == "chromosome"] <- "CHR"
names(results)[names(results) == "position"] <- "BP"
names(results)[names(results) == "rs_number"] <- "SNP"
names(results)[names(results) == "NMISS"] <- "NMISS"

if (!("CHR" %in% names(results)) || !("BP" %in% names(results))) {
  print("Chromosome and position not present in results dataframe. Pulling in now...")
  snps <- fread(snp_info, header = FALSE, stringsAsFactors = FALSE)
  snps <- snps[, c(1, 2, 4)]
  names(snps) <- c("CHR", "SNP", "BP")
  snps <- snps[snps$SNP %in% results$SNP, ]
  results <- merge(results, snps, by = "SNP", all.x = TRUE)
}

# Remove results with NA P values
if (sum(is.na(results$P)) > 0) {
  print(paste(sum(is.na(results$P)), "results with NA p values. Removing..."))
  results <- results[!is.na(results$P), ]
}

# Subset for viewing SNPs with P < 1e-5
significance <- subset(results, results$P < 1e-5)
write.table(significance, paste0(filename, "_significant"), quote = FALSE, row.names = FALSE)

# Calculate lambda
#chisq <- qchisq(1 - results$P, 1)
#lam <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

#print("Making QQ plot")
# QQ plot
#png(paste0(output_dir, basename(filename), ".qq.png"), type = "cairo")
#qq(results$P, sub = paste0("lambda = ", round(lam, 4)))
#dev.off()

# Remove results without CHR/BP info
if (sum(is.na(results$CHR)) > 0) {
  print(paste(sum(is.na(results$CHR)), "results without CHR/BP info. Removing to create the Manhattan plot..."))
  results <- results[!is.na(results$CHR), ]
}

print("Making Manhattan plot")
# Get results for Manhattan plot, including SNP column
man_results <- results[, c("CHR", "BP", "SNP", "P")]

# Define the SNPs to highlight
highlight_snps <- c("rs73072750", "rs28609108", "rs2959641")

# High-resolution TIFF plot for publication
tiff(paste0(output_dir, basename(filename), "_GrayScale.manhattan.tiff"), type = "cairo", 
     width = 3000, height = 1500, res = 300, compression = "lzw")

#manhattan(
 # man_results,
 # col = c("#000000", "#00003c", "#000079", "#0000b6", "#0000f2", "#1000fa", "#2600f3", "#3b00ed",
 #         "#5000e7", "#6c00c8", "#8e0098", "#af0067", "#d00037", "#e70f18", "#ed3612",
 #         "#f35e0b", "#fa8505", "#ffa900", "#ffbe00", "#ffd300", "#ffe900", "#ffff00"),
 # ylim = c(0, 10),
 # highlight = highlight_snps,
#  highlight.col = "red",
#  highlight.cex = 1.5
#)

manhattan(
  man_results,
  col = rep(c("black", "gray"), length.out = length(unique(man_results$CHR))),
  ylim = c(0, 10),
  highlight = highlight_snps,
  highlight.col = "red",
  highlight.cex = 1.5
)

dev.off()
print("Done with Manhattan plot!")