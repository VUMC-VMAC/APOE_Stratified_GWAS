args <- commandArgs(TRUE)
filename <- args[1] #include path with file stem
snp_info <- args[2] #bim file if these are meta-analysis results

library(qqman)
library(data.table)

# Define output directory
output_dir <- "/data/h_vmac/contra2/01_APOE4_Strat/03_GWAMA/plots/"

#read in results
results <- fread(filename, header=TRUE, stringsAsFactors = F)
names(results)[names(results) == "p-value"] <- "P"
names(results)[names(results) == "chromosome"] <- "CHR"
names(results)[names(results) == "position"] <- "BP"
names(results)[names(results) == "rs_number"] <- "SNP"
names(results)[names(results) == "NMISS"] <- "NMISS"

if(!("CHR" %in% names(results)) || !("BP" %in% names(results))){
  print("Chromosome and position not present in results dataframe. Pulling in now...")
  snps <- fread(snp_info, header = F, stringsAsFactors = F)
  snps <- snps[,c(1,2,4)]
  names(snps) <- c("CHR", "SNP", "BP")
  snps <- snps[snps$SNP %in% results$SNP,]
  results <- merge(results, snps, by = "SNP", all.x = T)
}

#remove results with NA
if(sum(is.na(results$P))>0){
  print(paste(sum(is.na(results$P)), "results with NA p values. Removing..."))
  results <- results[!is.na(results$P),]
}

##subset for viewing SNPs with p<5e-8
significance<-subset(results,results$P<0.00001)
write.table(significance,paste0(filename,"_significant"),quote=FALSE,row.names=FALSE)

#calculate lambda
chisq <- qchisq(1-results$P,1)
lam <- median(chisq,na.rm=TRUE)/qchisq(0.5,1)

print("Making qq plot")
#qqplot
png(paste0(output_dir, basename(filename), ".qq.png"), type="cairo")
qq(results$P, sub=paste0("lambda=", lam))
dev.off()

#remove results with NA
if(sum(is.na(results$CHR))>0){
  print(paste(sum(is.na(results$CHR)), "results without CHR/BP info. Removing to create the manhattan plot..."))
  results <- results[!is.na(results$CHR),]
}

print("Making manhattan plot")
#get results for manhattan plot
man_results <- results[,c("CHR","BP","P")]

#manhattan plot
png(paste0(output_dir, basename(filename), ".manhattan.png"), type="cairo", width = 960, height = 480)
manhattan(man_results,col=c("#000000","#00003c","#000079","#0000b6","#0000f2","#1000fa","#2600f3","#3b00ed",
"#5000e7","#6c00c8","#8e0098","#af0067","#d00037","#e70f18","#ed3612",
"#f35e0b","#fa8505","#ffa900","#ffbe00","#ffd300","#ffe900","#ffff00"),ylim=c(0,10))
dev.off()

print("Done with man_qq_plot!")