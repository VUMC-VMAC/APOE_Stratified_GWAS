
# module load GCC/11.3.0  OpenMPI/4.1.4
# module load R/4.2.1

# Directory containing the .GWAMA files
gwama_dir <- "/data/h_vmac/contra2/01_APOE4_Strat/02_GWAS/GWAMA_conversion"

library(data.table)
library(forestplot)

# List all .GWAMA files in the directory
gwama_files <- list.files(path = gwama_dir, pattern = "EXF.*APOE.*GWAMA", full.names = TRUE)

# Initialize an empty data frame to store results
rs2959641_data <- data.frame(cohort = character(), beta = numeric(), se = numeric(), n = integer(), group = character(), stringsAsFactors = FALSE)

# Loop over each file
for (file in gwama_files) {
  # Print a message indicating the file being read
  print(paste("Reading file:", file))
  
  # Use fread for faster reading
  gwama_data <- fread(file)
  
  # Extract the cohort name from the file path
  cohort_name <- tools::file_path_sans_ext(basename(file))
  
  # Modify cohort_name to shorten it, keeping only the first part
  cohort_name_short <- sub("^([^_]+_[^_]+).*", "\\1", cohort_name)
  
  # Determine if cohort is APOE4neg or APOE4pos based on file name
  group <- ifelse(grepl("APOE4neg", cohort_name), "APOE4neg", "APOE4pos")
  
  # Filter for rs2959641 and check if the SNP exists in the data
  rs2959641_row <- gwama_data[gwama_data$MARKERNAME == "rs2959641", ]
  
  if (nrow(rs2959641_row) > 0) {
    # Get beta, SE, and N for rs2959641
    beta <- rs2959641_row$BETA
    se <- rs2959641_row$SE
    n <- rs2959641_row$N
    
    # Append to rs2959641_data
    rs2959641_data <- rbind(rs2959641_data, data.frame(cohort = cohort_name_short, beta = beta, se = se, n = n, group = group))
  }
  
  # Print a message indicating completion of the file processing
  print(paste("Done with file:", file))
}


# View the extracted data
print(rs2959641_data)

# NHW with APOE4neg
rs2959641_NHW_neg <- subset(rs2959641_data, group == "APOE4neg" & grepl("_NHW", cohort))

# NHW with APOE4pos
rs2959641_NHW_pos <- subset(rs2959641_data, group == "APOE4pos" & grepl("_NHW", cohort))

# NHB with APOE4neg
rs2959641_NHB_neg <- subset(rs2959641_data, group == "APOE4neg" & grepl("_NHB", cohort))

# NHB with APOE4pos
rs2959641_NHB_pos <- subset(rs2959641_data, group == "APOE4pos" & grepl("_NHB", cohort))

#rs2959641_neg$weight <- 1 / (rs2959641_neg$se^2)
#rs2959641_pos$weight <- 1 / (rs2959641_pos$se^2)

########### RUNNING META ANALYSIS ###########
library(metafor)

meta_NHW_neg <- rma(yi = rs2959641_NHW_neg$beta, sei = rs2959641_NHW_neg$se)
meta_NHW_pos <- rma(yi = rs2959641_NHW_pos$beta, sei = rs2959641_NHW_pos$se)

meta_NHB_neg <- rma(yi = rs2959641_NHB_neg$beta, sei = rs2959641_NHB_neg$se)
meta_NHB_pos <- rma(yi = rs2959641_NHB_pos$beta, sei = rs2959641_NHB_pos$se)


# Pull coefficients
Beta_neg <- c(meta_NHW_neg$b, meta_NHB_neg$b)
SE_neg <- c(meta_NHW_neg$se, meta_NHB_neg$se)

Beta_pos <- c(meta_NHW_pos$b, meta_NHB_pos$b)
SE_pos <- c(meta_NHW_pos$se, meta_NHB_pos$se)

##Run cross-ancestry meta-analyses 
meta_neg <- rma(yi=SE_neg,sei=SE_neg,method="FE")
meta_pos <- rma(yi=Beta_pos,sei=SE_pos,method="FE")

#Major surgery

# Create data frame for APOE4 positive group (APOE4pos)
to_plot_APOE4pos <- as.data.frame(cbind(
  Beta = unlist(c(
    meta_NHW_pos$b[1], 
    subset(rs2959641_NHW_pos, cohort == "ACT_NHW")$beta, 
    subset(rs2959641_NHW_pos, cohort == "ADNI_NHW")$beta, 
    subset(rs2959641_NHW_pos, cohort == "BIOCARD_NHW")$beta, 
    subset(rs2959641_NHW_pos, cohort == "BLSA_NHW")$beta, 
    subset(rs2959641_NHW_pos, cohort == "NACC_NHW")$beta, 
    subset(rs2959641_NHW_pos, cohort == "ROSMAP_NHW")$beta, 
    subset(rs2959641_NHW_pos, cohort == "WASHU_NHW")$beta, 
    subset(rs2959641_NHW_pos, cohort == "WRAP_NHW")$beta, 
    meta_NHB_pos$b[1],
    subset(rs2959641_NHB_pos, cohort == "BLSA_NHB")$beta, 
    subset(rs2959641_NHB_pos, cohort == "NACC_NHB")$beta, 
    subset(rs2959641_NHB_pos, cohort == "ROSMAP_NHB")$beta, 
    meta_pos$b[1]
  )),
  
  lower = unlist(c(
    meta_NHW_pos$ci.lb[1], 
    subset(rs2959641_NHW_pos, cohort == "ACT_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "ACT_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "ADNI_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "ADNI_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "BIOCARD_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "BIOCARD_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "BLSA_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "BLSA_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "NACC_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "NACC_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "ROSMAP_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "ROSMAP_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "WASHU_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "WASHU_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "WRAP_NHW")$beta - 1.96 * subset(rs2959641_NHW_pos, cohort == "WRAP_NHW")$se,
    meta_NHB_pos$ci.lb,
    subset(rs2959641_NHB_pos, cohort == "BLSA_NHB")$beta - 1.96 * subset(rs2959641_NHB_pos, cohort == "BLSA_NHB")$se,
    subset(rs2959641_NHB_pos, cohort == "NACC_NHB")$beta - 1.96 * subset(rs2959641_NHB_pos, cohort == "NACC_NHB")$se,
    subset(rs2959641_NHB_pos, cohort == "ROSMAP_NHB")$beta - 1.96 * subset(rs2959641_NHB_pos, cohort == "ROSMAP_NHB")$se,
    meta_pos$ci.lb[1]
  )),
  
  upper = unlist(c(
    meta_NHW_pos$ci.ub[1], 
    subset(rs2959641_NHW_pos, cohort == "ACT_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "ACT_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "ADNI_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "ADNI_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "BIOCARD_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "BIOCARD_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "BLSA_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "BLSA_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "NACC_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "NACC_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "ROSMAP_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "ROSMAP_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "WASHU_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "WASHU_NHW")$se,
    subset(rs2959641_NHW_pos, cohort == "WRAP_NHW")$beta + 1.96 * subset(rs2959641_NHW_pos, cohort == "WRAP_NHW")$se,
    meta_NHW_pos$ci.ub,
    subset(rs2959641_NHB_pos, cohort == "BLSA_NHB")$beta + 1.96 * subset(rs2959641_NHB_pos, cohort == "BLSA_NHB")$se,
    subset(rs2959641_NHB_pos, cohort == "NACC_NHB")$beta + 1.96 * subset(rs2959641_NHB_pos, cohort == "NACC_NHB")$se,
    subset(rs2959641_NHB_pos, cohort == "ROSMAP_NHB")$beta + 1.96 * subset(rs2959641_NHB_pos, cohort == "ROSMAP_NHB")$se,
    meta_pos$ci.ub[1]
  ))
))

# Create data frame for APOE4 negitive group (APOE4neg)
to_plot_APOE4neg <- as.data.frame(cbind(
  Beta = unlist(c(
    meta_NHW_neg$b[1], 
    subset(rs2959641_NHW_neg, cohort == "ACT_NHW")$beta, 
    subset(rs2959641_NHW_neg, cohort == "ADNI_NHW")$beta, 
    subset(rs2959641_NHW_neg, cohort == "BIOCARD_NHW")$beta, 
    subset(rs2959641_NHW_neg, cohort == "BLSA_NHW")$beta, 
    subset(rs2959641_NHW_neg, cohort == "NACC_NHW")$beta, 
    subset(rs2959641_NHW_neg, cohort == "ROSMAP_NHW")$beta, 
    subset(rs2959641_NHW_neg, cohort == "WASHU_NHW")$beta, 
    subset(rs2959641_NHW_neg, cohort == "WRAP_NHW")$beta, 
    meta_NHB_neg$b[1],
    subset(rs2959641_NHB_neg, cohort == "BLSA_NHB")$beta, 
    subset(rs2959641_NHB_neg, cohort == "NACC_NHB")$beta, 
    subset(rs2959641_NHB_neg, cohort == "ROSMAP_NHB")$beta, 
    meta_neg$b[1]
  )),
  
  lower = unlist(c(
    meta_NHW_neg$ci.lb[1], 
    subset(rs2959641_NHW_neg, cohort == "ACT_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "ACT_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "ADNI_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "ADNI_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "BIOCARD_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "BIOCARD_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "BLSA_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "BLSA_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "NACC_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "NACC_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "ROSMAP_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "ROSMAP_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "WASHU_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "WASHU_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "WRAP_NHW")$beta - 1.96 * subset(rs2959641_NHW_neg, cohort == "WRAP_NHW")$se,
    meta_NHB_neg$ci.lb,
    subset(rs2959641_NHB_neg, cohort == "BLSA_NHB")$beta - 1.96 * subset(rs2959641_NHB_neg, cohort == "BLSA_NHB")$se,
    subset(rs2959641_NHB_neg, cohort == "NACC_NHB")$beta - 1.96 * subset(rs2959641_NHB_neg, cohort == "NACC_NHB")$se,
    subset(rs2959641_NHB_neg, cohort == "ROSMAP_NHB")$beta - 1.96 * subset(rs2959641_NHB_neg, cohort == "ROSMAP_NHB")$se,
    meta_neg$ci.lb[1]
  )),
  
  upper = unlist(c(
    meta_NHW_neg$ci.ub[1], 
    subset(rs2959641_NHW_neg, cohort == "ACT_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "ACT_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "ADNI_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "ADNI_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "BIOCARD_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "BIOCARD_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "BLSA_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "BLSA_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "NACC_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "NACC_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "ROSMAP_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "ROSMAP_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "WASHU_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "WASHU_NHW")$se,
    subset(rs2959641_NHW_neg, cohort == "WRAP_NHW")$beta + 1.96 * subset(rs2959641_NHW_neg, cohort == "WRAP_NHW")$se,
    meta_NHW_neg$ci.ub,
    subset(rs2959641_NHB_neg, cohort == "BLSA_NHB")$beta + 1.96 * subset(rs2959641_NHB_neg, cohort == "BLSA_NHB")$se,
    subset(rs2959641_NHB_neg, cohort == "NACC_NHB")$beta + 1.96 * subset(rs2959641_NHB_neg, cohort == "NACC_NHB")$se,
    subset(rs2959641_NHB_neg, cohort == "ROSMAP_NHB")$beta + 1.96 * subset(rs2959641_NHB_neg, cohort == "ROSMAP_NHB")$se,
    meta_neg$ci.ub[1]
  ))
))


#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@#@@@@@@@@@@@@
#      GCC/11.3.0  OpenMPI/4.1.4
# module load R/4.2.1
tabletext <- list(c("Non-Hispanic White", "         ACT","         ADNI","         BIOCARD","         BLSA","         NACC", "         ROS/MAP/MARS","         WASHU","         WRAP",  "Non-Hispanic Black", "         BLSA", "         NACC", "         ROS/MAP/MARS", "Cross-Ancestry"))

dir <- "/data/h_vmac/contra2/01_APOE4_Strat/07_Other/Forest_Plots"

tiff(paste0(dir,"/Forest_plot_rs2959641.tiff"), width=8, height=7, units="in", res=300, type="cairo")
forestplot(tabletext, 
           mean = cbind(to_plot_APOE4pos$Beta, to_plot_APOE4neg$Beta), 
           is.summary = c(TRUE, rep(FALSE, 8), TRUE, rep(FALSE, 3), TRUE),
           col = fpColors(box = c("mediumseagreen", "mediumpurple"), line = c("mediumseagreen", "mediumpurple"),
                          summary = c("mediumseagreen", "mediumpurple")),
           lower = cbind(to_plot_APOE4pos$lower, to_plot_APOE4neg$lower),
           upper = cbind(to_plot_APOE4pos$upper, to_plot_APOE4neg$upper), 
           xlab = "Beta [95% CI]", ci.vertices = TRUE,
           xticks = seq(-0.3, 0.3, by = 0.1),  # Adjust xticks to match new x-axis range
           xlim = c(-0.3, 0.3),                # Set x-axis limits           ci.vertices.height = 0.07, boxsize = 0.1, 
           legend = c(expression(italic("APOE") * "-ε4 Carriers"), "Non-carriers"), 
           title = "rs2959641 on Baseline Executive Function")

dev.off()