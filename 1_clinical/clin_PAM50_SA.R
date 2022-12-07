# Script: Survival analyses based on second best pam50

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_PAM50/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/1_clinical/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/1_clinical/src/pam50_functions.R")
library(ggplot2)
library(ggfortify)
library(survival)
library(tidyverse)
library(survminer)
library(grid)

# list to store plots
plot.list <- list()

#######################################################################
# Load data
#######################################################################

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData") %>% 
  filter(EvalGroup_ERpHER2nLNn50 == "ERpHER2nLNn_Endo_50") %>% 
  filter(majorityNearestClass != majoritySecondBestClass)

pam50.subtypes <- unique(anno$NCN_PAM50)
outcome.measures <- c("OS","DRFI")

#test
pam50.type <- "LumA"
data <- anno[which(anno$NCN_PAM50 == pam50.type),]
OM <- "OS"
source("scripts/1_clinical/src/pam50_functions.R")


# create KM plots comparing OS and DFRI per subtype per second best matching subtype
for (j in 1:length(pam50.subtypes)) {
  # select data
  pam50.type <- pam50.subtypes[j]
  pam50.data <- anno[which(anno$NCN_PAM50 == pam50.type),]
  # do SA based on majoritySecondBestClass
  for (k in 1:length(outcome.measures)) {
    OM <- outcome.measures[k]
    OM.bin <- paste(OM, "bin",sep="")    
    # surv object
    data.surv <- Surv(pam50.data[[OM]], pam50.data[[OM.bin]]) 
    # fit
    fit <- survminer::surv_fit(data.surv~majoritySecondBestClass, data=pam50.data, conf.type="log-log")
    # plot
    plot.list <- append(plot.list,list(
      KM.plot(pam50.data,fit,OM,
                   colors = c("#d334eb", "#2176d5", "#34c6eb"),
                   title = paste(OM, " in NCN_PAM50: ", pam50.type, sep=""),
                   group.variable = "majoritySecondBestClass",
                   legend.labels = c("HER2E","LUMB","NORMAL"))))
    
  }
  
}


    
    plot <- ggsurvplot(
        fit,
        censor.size = 8,
        censor.shape = "|",
        size = 5,
        risk.table = FALSE,       
        pval = TRUE,
        pval.size = 8,
        pval.coord = c(0,0.1),
        conf.int = FALSE,         
        xlim = c(0,max(clin.rel4$OS_rel4[is.finite(clin.rel4$OS_rel4)])),         
        xlab = "OS days",
        ylab = "OS event probability", # ggf just label as "event probability"
        ylim = c(0,1),
        palette = c("#d334eb", "#2176d5", "#34c6eb"), 
        legend = c(0.9,0.96),
        ggtheme = theme(legend.title = element_text(size=25), #20
                        legend.key.size = unit(0.5,"cm"), 
                        legend.text = element_text(size = 25), #20
                        axis.text.x = element_text(size = 25), #20
                        axis.title.x = element_text(size = 30), #25
                        axis.text.y = element_text(size = 25), #20
                        axis.title.y = element_text(size = 30),
                        plot.title = element_text(size=30)),
        title= "Overall survival in FGFR4 expression groups",
        legend.title = "FGFR4",
        legend.labs = c(paste("High"," (",table(clin.rel4[!is.na(clin.rel4$OS_rel4),]$FGFR4_expr)[1],")",sep = ""),
                        paste("Low"," (",table(clin.rel4[!is.na(clin.rel4$OS_rel4),]$FGFR4_expr)[2],")",sep = ""),
                        paste("Medium"," (",table(clin.rel4[!is.na(clin.rel4$OS_rel4),]$FGFR4_expr)[3],")",sep = "")),
        break.x.by = 500, # break X axis in time intervals of x (what is nicest here? maybe 365)
        break.y.by = 0.1)
    
    print(plot)
    
# RFI
# surv object
data.surv <- Surv(clin.rel4$RFI_rel4, clin.rel4$RFIbin_rel4) 
    
# fit
fit <- survminer::surv_fit(data.surv~FGFR4_expr, data=clin.rel4, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
      fit,
      censor.size = 8,
      censor.shape = "|",
      size = 5,
      risk.table = FALSE,       
      pval = TRUE,
      pval.size = 8,
      pval.coord = c(0,0.1),
      conf.int = FALSE,         
      xlim = c(0,max(clin.rel4$RFI_rel4[is.finite(clin.rel4$RFI_rel4)])),         
      xlab = "RFI days",
      ylab = "RFI event probability", # ggf just label as "event probability"
      ylim = c(0,1),
      palette = c("#d334eb", "#2176d5", "#34c6eb"), 
      legend = c(0.9,0.96),
      ggtheme = theme(legend.title = element_text(size=25), #20
                      legend.key.size = unit(0.5,"cm"), 
                      legend.text = element_text(size = 25), #20
                      axis.text.x = element_text(size = 25), #20
                      axis.title.x = element_text(size = 30), #25
                      axis.text.y = element_text(size = 25), #20
                      axis.title.y = element_text(size = 30),
                      plot.title = element_text(size=30)),
      title= "RFI in FGFR4 expression groups",
      legend.title = "FGFR4",
      legend.labs = c(paste("High"," (",table(clin.rel4[!is.na(clin.rel4$RFI_rel4),]$FGFR4_expr)[1],")",sep = ""),
                      paste("Low"," (",table(clin.rel4[!is.na(clin.rel4$RFI_rel4),]$FGFR4_expr)[2],")",sep = ""),
                      paste("Medium"," (",table(clin.rel4[!is.na(clin.rel4$RFI_rel4),]$FGFR4_expr)[3],")",sep = "")),
      break.x.by = 500, # break X axis in time intervals of x (what is nicest here? maybe 365)
      break.y.by = 0.1)
    
    print(plot)
    

dev.off()

