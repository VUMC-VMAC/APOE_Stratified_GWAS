#Build covariate/phenotype files for global cognition for multi-domain analyses 
library(dplyr)

Main_directory<-"/data/h_vmac/contra2/01_APOE4_Strat/01_CovPheno/0_CognitiveBuild/"

#Setting directory names
File_directory <- paste0(Main_directory, "GWAS_Files/")
NHB_IDs_directory <- paste0(Main_directory, "NHB_IDs/")


#reading cognitive data
cogn_data <- readRDS(paste0(Main_directory, "Cognitive_Master_Build_June_2023_v2.RDS"))


##########################################################################################
###################Create phenotype and covariate files for MEM/EXF/LAN###################
##########################################################################################
cog_variables <- c("MEM", "memslopes", "EXF", "exfslopes", "LAN", "lanslopes")
cog_dictionary <- c("memslopes" = "MEM", "MEM" = "MEM","exfslopes" = "EXF", "EXF" = "EXF", "lanslopes" = "LAN", "LAN" = "LAN")
cohorts <- c("ACT", "ADNI", "BIOCARD", "BLSA", "NACC", "ROSMAP", "WASHU", "WRAP") 
pc_dictionary <- c("NHW"="NHW","NHB"="NHB")
relateds_dictionary <- c("NHW"="NHW","NHB"="AllRaces")

for (i in 1:length(cog_variables)){
  
  for (j in 1:length(cohorts)){
    
    ###Setup variables.
    cohort <- cohorts[j]
    var <- cog_variables[i]
    interval_var <- paste0("interval_", cog_dictionary[var])
    Age_var <- paste0("Age_bl_", cog_dictionary[var])
    dx_var <- paste0("dx_bl_", cog_dictionary[var])
    num_visits_var <- paste0("num_visits_", cog_dictionary[var])
    
    ##################################
    #####Save all necessary files#####
    ##################################
    race_categories <- c("NHW","NHB")
    
    for (k in 1:length(race_categories)){
      
      print(k)
      ###Get baseline data point for cognitive measure. 
      data_race <- cogn_data[(cogn_data$Age == cogn_data[Age_var]),] %>% 
        filter(Study == cohort) %>% 
        filter(race == race_categories[k])
      
      
      ###Remove all relateds.
      {
        participants_to_include <- read.table(paste0(Main_directory, "NonRelated_Participants/", relateds_dictionary[k], "/", cohort, "_", relateds_dictionary[k], "_imputed_final_no_relateds.txt"), header=TRUE) 
        
        ###Create unique FID variable.
        data_race$FID_races <- paste0(cohort, "_", pc_dictionary[k], "_", data_race$FID)
        
        ###Remove participants not in participants_to_include dataframe.
        print(paste0(nrow(data_race), " participants for ", cohort, " before ", race_categories[k], " relatedness removal!"))
        data_race <- data_race[(data_race$FID %in% participants_to_include$FID),]
        print(paste0(nrow(data_race), " participants for ", cohort, " after ", race_categories[k], " relatedness removal!"))
        }
      
      ##################################
      #######Get and save NHB IDs#######
      ##################################
      if(race_categories[k] == "NHB"){
        NHB_IDs <- data.frame(FID = data_race$FID, IID = data_race$IID)
        write.table(NHB_IDs, file=paste0(NHB_IDs_directory, cohort, "_", race_categories[k], "_","IDs.txt"), row.names=FALSE, quote=FALSE)
      }
      
      
      #########################################
      #####Create and save phenotype files#####
      #########################################
      
      ###Setup phenotype files 
      data_phenotype_with_comorbid_participants <- data_race %>% 
        filter(!is.na(!!(sym(var)))) %>% 
        filter(!is.na(IID)) %>% 
        filter(!is.na(!!as.name(paste0(pc_dictionary[k], "_PC1")))) 

      ###Save 50+ phenotype files
      if(nrow(data_phenotype_with_comorbid_participants %>% filter(Age>=50))>1){write.table(data_phenotype_with_comorbid_participants %>% filter(Age>=50) %>% dplyr::select(c("FID", "IID", !!(sym(var)))), file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus_ALL", ".txt"), row.names=FALSE, quote=FALSE)}
      #########################################
      #####Create and save covariate files#####
      #########################################
      
      ###PC variable names
      PC1 <- paste0(race_categories[k], "_PC1")
      PC2 <- paste0(race_categories[k], "_PC2")
      PC3 <- paste0(race_categories[k], "_PC3")
      PC4 <- paste0(race_categories[k], "_PC4")
      PC5 <- paste0(race_categories[k], "_PC5")
      
      ###Setup cross-sectional covariate files. 
      covariate_variable_names <-c("FID", "IID", "Age", "sex", "education", "race", "dx",
                                   "apoe4count", "apoe2count", "apoe4pos", "apoe2pos",
                                   PC1, PC2, PC3, PC4, PC5, var, "Comorbidities")
      
     data_cov_with_comorbid_participants_50plus <- data_race[covariate_variable_names] %>% filter(!is.na(PC1)) %>%  filter(!is.na(!!(sym(var)))) %>% filter(!is.na(IID))
     
      ###Save 50+ covariate files
      ###ALL
      if(nrow(data_cov_with_comorbid_participants_50plus)>1){write.table(data_cov_with_comorbid_participants_50plus, file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus", "_", "MainEffects_covar_ALL", ".txt"), row.names=FALSE, quote=FALSE)}
      print(paste0(cohorts[j], ":", race_categories[k], ":", var, ":", "MainEffects:", "WithComorbidities", ":50plus", ": ", nrow(data_cov_with_comorbid_participants_50plus)))

      ###APOE4neg
      if(nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==0))>1){write.table(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==0),file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus", "_", "MainEffects_covar_APOE4neg", ".txt"),row.names=FALSE,quote=FALSE)}
      print(paste0(cohorts[j], ":", race_categories[k], ":", var, ":", "APOE4neg:Maineffects:", "WithComorbidities", ":50plus", ": " ,nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==0))))

      ###APOE4pos
      if(nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==1))>1){write.table(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==1),file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus", "_", "MainEffects_covar_APOE4pos", ".txt"),row.names=FALSE,quote=FALSE)}
      print(paste0(cohorts[j], ":", race_categories[k], ":", var, ":", "APOE4pos:MainEffects:", "WithComorbidities", ":50plus", ": ", nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==1))))

    }
  } 
}


#Doing the same thing but for just A4
##########################################################################################
###################Create phenotype and covariate files for MEM/EXF/LAN###################
##########################################################################################
cog_variables <- c("MEM")
cog_dictionary <- c("MEM" = "MEM")
cohorts <- c("A4") 
pc_dictionary <- c("NHW"="NHW","NHB"="NHB")
relateds_dictionary <- c("NHW"="NHW","NHB"="AllRaces")

for (i in 1:length(cog_variables)){
  
  for (j in 1:length(cohorts)){
    
    ###Setup variables.
    cohort <- cohorts[j]
    var <- cog_variables[i]
    interval_var <- paste0("interval_", cog_dictionary[var])
    Age_var <- paste0("Age_bl_", cog_dictionary[var])
    dx_var <- paste0("dx_bl_", cog_dictionary[var])
    num_visits_var <- paste0("num_visits_", cog_dictionary[var])
    
    ##################################
    #####Save all necessary files#####
    ##################################
    race_categories <- c("NHW","NHB")
    
    for (k in 1:length(race_categories)){
      
      print(k)
      ###Get baseline data point for cognitive measure. 
      data_race <- cogn_data[(cogn_data$Age == cogn_data[Age_var]),] %>% 
        filter(Study == cohort) %>% 
        filter(race == race_categories[k])
      
      
      ###Remove all relateds.
      {
        participants_to_include <- read.table(paste0(Main_directory, "NonRelated_Participants/", relateds_dictionary[k], "/", cohort, "_", relateds_dictionary[k], "_imputed_final_no_relateds.txt"), header=TRUE) 
        
        ###Create unique FID variable.
        data_race$FID_races <- paste0(cohort, "_", pc_dictionary[k], "_", data_race$FID)
        
        ###Remove participants not in participants_to_include dataframe.
        print(paste0(nrow(data_race), " participants for ", cohort, " before ", race_categories[k], " relatedness removal!"))
        data_race <- data_race[(data_race$FID %in% participants_to_include$FID),]
        print(paste0(nrow(data_race), " participants for ", cohort, " after ", race_categories[k], " relatedness removal!"))
        }
      
      ##################################
      #######Get and save NHB IDs#######
      ##################################
      if(race_categories[k] == "NHB"){
        NHB_IDs <- data.frame(FID = data_race$FID, IID = data_race$IID)
        write.table(NHB_IDs, file=paste0(NHB_IDs_directory, cohort, "_", race_categories[k], "_","IDs.txt"), row.names=FALSE, quote=FALSE)
      }
      
      
      #########################################
      #####Create and save phenotype files#####
      #########################################
      
      ###Setup phenotype files 
      data_phenotype_with_comorbid_participants <- data_race %>% 
        filter(!is.na(!!(sym(var)))) %>% 
        filter(!is.na(IID)) %>% 
        filter(!is.na(!!as.name(paste0(pc_dictionary[k], "_PC1")))) 
      
      ###Save 50+ phenotype files
      if(nrow(data_phenotype_with_comorbid_participants %>% filter(Age>=50))>1){write.table(data_phenotype_with_comorbid_participants %>% filter(Age>=50) %>% dplyr::select(c("FID", "IID", !!(sym(var)))), file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus_ALL", ".txt"), row.names=FALSE, quote=FALSE)}
      #########################################
      #####Create and save covariate files#####
      #########################################
      
      ###PC variable names
      PC1 <- paste0(race_categories[k], "_PC1")
      PC2 <- paste0(race_categories[k], "_PC2")
      PC3 <- paste0(race_categories[k], "_PC3")
      PC4 <- paste0(race_categories[k], "_PC4")
      PC5 <- paste0(race_categories[k], "_PC5")
      
      ###Setup cross-sectional covariate files. 
      covariate_variable_names <-c("FID", "IID", "Age", "sex", "education", "race", "dx",
                                   "apoe4count", "apoe2count", "apoe4pos", "apoe2pos",
                                   PC1, PC2, PC3, PC4, PC5, var, "Comorbidities")
      
      data_cov_with_comorbid_participants_50plus <- data_race[covariate_variable_names] %>% filter(!is.na(PC1)) %>%  filter(!is.na(!!(sym(var)))) %>% filter(!is.na(IID))
      
      ###Save 50+ covariate files
      ###ALL
      if(nrow(data_cov_with_comorbid_participants_50plus)>1){write.table(data_cov_with_comorbid_participants_50plus, file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus", "_", "MainEffects_covar_ALL", ".txt"), row.names=FALSE, quote=FALSE)}
      print(paste0(cohorts[j], ":", race_categories[k], ":", var, ":", "MainEffects:", "WithComorbidities", ":50plus", ": ", nrow(data_cov_with_comorbid_participants_50plus)))
      
      ###APOE4neg
      if(nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==0))>1){write.table(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==0),file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus", "_", "MainEffects_covar_APOE4neg", ".txt"),row.names=FALSE,quote=FALSE)}
      print(paste0(cohorts[j], ":", race_categories[k], ":", var, ":", "APOE4neg:Maineffects:", "WithComorbidities", ":50plus", ": " ,nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==0))))
      
      ###APOE4pos
      if(nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==1))>1){write.table(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==1),file=paste0(File_directory, cohort, "_", race_categories[k], "_", var, "_", "WithComorbidities_50plus", "_", "MainEffects_covar_APOE4pos", ".txt"),row.names=FALSE,quote=FALSE)}
      print(paste0(cohorts[j], ":", race_categories[k], ":", var, ":", "APOE4pos:MainEffects:", "WithComorbidities", ":50plus", ": ", nrow(data_cov_with_comorbid_participants_50plus %>% filter(apoe4pos==1))))
      
    }
  } 
}






































