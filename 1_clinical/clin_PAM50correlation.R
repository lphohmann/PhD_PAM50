# Script: Checking PAM50 correlation

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_PAM50/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- "data/SCANB/1_clinical/processed/"
dir.create(data.path)

#packages
source("./scripts/1_clinical/src/pam50_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
library(openxlsx)
library(readxl)
library(ggfortify)
library(survival)
library(survminer)
library(grid)

# list to store plots
plot.list <- list()

#######################################################################
# Load data
#######################################################################

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData")

pam50.subtypes <- unique(anno$NCN_PAM50)

clinical.groups <- setNames(as.data.frame(rbind(
  All = c(NA,NA), 
  TNBC = c("Negative","Negative"), 
  ERpHER2p = c("Positive","Positive"), 
  ERnHER2p = c("Negative","Positive"), 
  ERpHER2n = c("Positive","Negative")), 
  stringsAsFactors = FALSE), 
  c("ER","HER2"))

#######################################################################
# PAM50 centroid correlation
#######################################################################

for (i in 1:nrow(clinical.groups)) { 
  # select subgroup data
  clin.group <- row.names(clinical.groups)[i] #i
  if (clin.group == "All") {
    sg.data <- anno
  } else {
    sg.data <- anno[which(
      anno$HER2 == clinical.groups[clin.group ,"HER2"] & # check format in anno
      anno$ER == clinical.groups[clin.group ,"ER"]),]
  }
  for (j in 1:length(pam50.subtypes)) {
    # select data 
    pam50.type <- pam50.subtypes[j]
    pam50.data <- sg.data[which(sg.data$NCN_PAM50 == pam50.type),]
    # do
    pam50.data.mod <- melt(pam50.data, measure.vars=colnames(pam50.data[grepl("mean",colnames(pam50.data))])) # only work if other columns dont caintain "mean" substring
    plot.list <- append(plot.list,list(
      ggplot(pam50.data.mod) + 
        geom_boxplot(aes(x=variable, y=value,fill=as.factor(variable)),alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("PAM50 Centroid") +
        ylab("Correlation") +
        ggtitle(paste(
          "PAM50 centroid correlation - clin: ",clin.group,"; pam50: ",pam50.type,sep="")) +
        theme(plot.title = element_text(size = 25),
              axis.text.x = element_text(size = 20),
              axis.title.x = element_text(size = 25),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 25),
              legend.position = "none") +
        scale_fill_manual(values=group.colors <- c(meanLumA = "#2176d5", meanLumB = "#34c6eb", meanHer2 ="#d334eb", meanBasal ="#c41b0e", meanNormal="#64c25d"))))
  }
}

#######################################################################
# Crosstable with 2nd best PAM50 match for all samples
#######################################################################

for (i in 1:nrow(clinical.groups)) { 
  # select subgroup data
  clin.group <- row.names(clinical.groups)[i] #i
  anno <- anno[which(anno$majorityNearestClass!=anno$majoritySecondBestClass),]
  if (clin.group == "All") {
    sg.data <- anno
  } else {
    sg.data <- anno[which(
      anno$HER2 == clinical.groups[clin.group ,"HER2"] & 
        anno$ER == clinical.groups[clin.group ,"ER"]),]
  }
  
  # get cm
  confusion_matrix <- as.data.frame(table(
    sg.data$majorityNearestClass,sg.data$majoritySecondBestClass,dnn=c("majorityNearestClass","majoritySecondBestClass")))
  
  # plot
  plot.list <- append(plot.list,list(
    ggplot(data = confusion_matrix,
           mapping = aes(x = majorityNearestClass,
                         y = majoritySecondBestClass)) +
      scale_x_discrete(limits = rev(levels(confusion_matrix$majorityNearestClass))) +
      geom_tile(aes(fill = Freq)) +
      geom_text(aes(label = sprintf("%1.0f", Freq)), size=10, vjust = 1) +
      scale_fill_gradient(low = "#fee6ce",
                          high = "#e6550d",
                          trans = "log") +
      ggtitle(paste("Second-best PAM50 match - clin: ",clin.group,sep="")) +
      xlab("majorityNearestClass") +
      ylab("majoritySecondBestClass") +
      theme(plot.title = element_text(size = 25),
            axis.text.x = element_text(size = 20),
            axis.title.x = element_text(size = 25),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 25), 
            legend.position = "none")))
}


#######################################################################
# Distinctiveness: Boxplot the difference between the best and second 
# best correlation matches
#######################################################################

for (i in 1:nrow(clinical.groups)) { 
  # select subgroup data
  clin.group <- row.names(clinical.groups)[i] 
  # filter out samples where best=2nd best and calc diff.
  anno <- anno[which(
    anno$majorityNearestClass!=anno$majoritySecondBestClass),] %>% 
    rowwise() %>% mutate(Diff = abs(get(paste("mean",majorityNearestClass,sep=""))) - abs(get(paste("mean",majoritySecondBestClass,sep=""))))
  
  if (clin.group == "All") {
    sg.data <- anno
  } else {
    sg.data <- anno[which(
      anno$HER2 == clinical.groups[clin.group ,"HER2"] & 
        anno$ER == clinical.groups[clin.group ,"ER"]),]
  }
  
  # plot
  plot.list <- append(plot.list,list(
    ggplot(sg.data) + 
    geom_boxplot(aes(x=majorityNearestClass, y=Diff,fill=as.factor(majorityNearestClass)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 Class") +
    ylab("Correlation diff (1st-2nd)") +
    ggtitle(paste("PAM50 class distinctiveness - clin: ",clin.group,sep="")) +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_fill_manual(values= c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb", Basal ="#c41b0e", Normal="#64c25d"))))
}

#######################################################################
# PAM50 Class specific distinctiveness: difference between the best and second best correlation matches
#######################################################################

for (i in 1:nrow(clinical.groups)) { 
  # select subgroup data
  clin.group <- row.names(clinical.groups)[i] 
  # filter out samples where best=2nd best and calc diff.
  anno <- anno[which(
    anno$majorityNearestClass!=anno$majoritySecondBestClass),] %>% 
    rowwise() %>% mutate(Diff = abs(get(paste("mean",majorityNearestClass,sep=""))) - abs(get(paste("mean",majoritySecondBestClass,sep=""))))
  
  if (clin.group == "All") {
    sg.data <- anno
  } else {
    sg.data <- anno[which(
      anno$HER2 == clinical.groups[clin.group ,"HER2"] & 
        anno$ER == clinical.groups[clin.group ,"ER"]),]
  }
  
  for (j in 1:length(pam50.subtypes)) {
    # select data 
    pam50.type <- pam50.subtypes[j]
    pam50.data <- sg.data[which(sg.data$NCN_PAM50 == pam50.type),]
    
    # plot
    plot.list <- append(plot.list,list(
      ggplot(pam50.data) + 
        geom_boxplot(aes(x=majoritySecondBestClass, y=Diff,fill=as.factor(majoritySecondBestClass)),alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("majoritySecondBestClass") +
        ylab("Correlation diff (1st-2nd)") +
        ggtitle(paste("PAM50 class-specific distinctiveness - clin: ",clin.group,"; pam50: ",pam50.type,sep="")) +
        theme(plot.title = element_text(size = 25),
              axis.text.x = element_text(size = 20),
              axis.title.x = element_text(size = 25),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 25),
              legend.position = "none") +
        scale_fill_manual(values= c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb", Basal ="#c41b0e", Normal="#64c25d"))))
  }
}

#######################################################################
# SA for second best pam50 class
#######################################################################

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData") %>% 
  filter(EvalGroup_ERpHER2nLNn50 == "ERpHER2nLNn_Endo_50") %>% 
  filter(majorityNearestClass != majoritySecondBestClass)

pam50.subtypes <- unique(anno$NCN_PAM50)
outcome.measures <- c("OS","DRFI")
colors <- c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb", Basal ="#c41b0e", Normal = "#64c25d")

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
              colors = unname(colors[names(table(pam50.data$majoritySecondBestClass))]),
              title = paste(OM, " in ERpHER2nLNn_Endo_50; NCN_PAM50: ", pam50.type, sep=""),
              group.variable = "majoritySecondBestClass",
              legend.labels = names(table(pam50.data$majoritySecondBestClass))
              ))) 
  }
  #colors[grepl(paste(names(table(pam50.data$majoritySecondBestClass)), collapse = '|'),names(colors))],
}

#plot
pdf(file = paste(output.path,"PAM50_analyses.pdf", sep=""), 
    onefile = TRUE, width = 21.0, height = 14.8) 

for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}

dev.off()


#######################################################################
# SA for disctinctiveness in LumA and LumB
#######################################################################

# see what is worth creating a function for

# distinctiveness based on median distinctiveness in these subgroups (e.g. LumA-Normal -> split based on median)

# LumA

# KM plot: LumA-Normal (distinct); LumA-Normal (non-distinct); LumA-LumB (distinct); LumA-LumB (non-distinct)

# Boxplots: Age; Tum Size (with these groups)

# LumB

# KM plot: LumB-Normal (distinct); LumB-Normal (non-distinct); LumB-LumA (distinct); LumB-LumA (non-distinct)

# Boxplots: Age; Tum Size (with these groups)






#######################################################################
# plot correlations
#######################################################################
# 
# # all samples
# plot.list <- append(plot.list,list(
#     ggplot(data_mod) + 
#     geom_point(aes(x=meanHer2, y=HER2Diff,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
#     xlab("HER2E correlation (mean)") +
#     ylab("Correlation diff (1st-2nd)") +
#     ggtitle("HER2E distinctiveness: Subtype centroid correlations (ERpHER2nHER2E)") +
#     theme(plot.title = element_text(size = 20),
#           axis.text.x = element_text(size = 20),
#           axis.title.x = element_text(size = 25),
#           axis.text.y = element_text(size = 20),
#           axis.title.y = element_text(size = 25),
#           legend.position = c(0.9,0.1),
#           legend.text=element_text(size = 15),
#           legend.title=element_text(size = 20)) +
#     scale_color_manual("2nd Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
#     scale_x_continuous(breaks=seq(0, 0.8, 0.1))))
# 
# # correlation between subtype and her2 centroid correlations
# # luma
# data_mod <- data %>% filter(majoritySecondBest == "LumA")
# 
# plot.list <- append(plot.list,list(
#     ggplot(data_mod) + 
#     geom_point(aes(x=meanHer2, y=meanLumA,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
#     xlab("HER2E correlation") +
#     ylab("LumA correlation") +
#     ggtitle("ERpHER2nHER2E: Luminal A and HER2E centroid correlations") +
#     theme(plot.title = element_text(size = 20),
#           axis.text.x = element_text(size = 20),
#           axis.title.x = element_text(size = 25),
#           axis.text.y = element_text(size = 20),
#           axis.title.y = element_text(size = 25),
#           legend.position = "none") +
#     scale_color_manual("Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
#     scale_x_continuous(limits = c(0.2, 0.8),breaks=seq(0, 0.8, 0.1)) +
#     scale_y_continuous(limits = c(0,0.6),breaks=seq(0, 0.6, 0.1))))
# 
# # lumb
# data_mod <- data %>% filter(majoritySecondBest == "LumB")
# 
# plot.list <- append(plot.list,list(
#     ggplot(data_mod) + 
#     geom_point(aes(x=meanHer2, y=meanLumB,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
#     xlab("HER2E correlation") +
#     ylab("LumB correlation") +
#     ggtitle("ERpHER2nHER2E: Luminal B and HER2E centroid correlations") +
#     theme(plot.title = element_text(size = 20),
#           axis.text.x = element_text(size = 20),
#           axis.title.x = element_text(size = 25),
#           axis.text.y = element_text(size = 20),
#           axis.title.y = element_text(size = 25),
#           legend.position = "none") +
#     scale_color_manual("Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
#     scale_x_continuous(limits = c(0.2, 0.8),breaks=seq(0, 0.8, 0.1)) +
#     scale_y_continuous(limits = c(0,0.6),breaks=seq(0, 0.6, 0.1))))
# 
# # basal
# data_mod <- data %>% filter(majoritySecondBest == "Basal")
# 
# plot.list <- append(plot.list,list(
#     ggplot(data_mod) + 
#     geom_point(aes(x=meanHer2, y=meanBasal,color=as.factor(majoritySecondBest)),alpha=0.7, size=5) +
#     xlab("HER2E correlation") +
#     ylab("Basal correlation") +
#     ggtitle("ERpHER2nHER2E: Basal and HER2E centroid correlations") +
#     theme(plot.title = element_text(size = 20),
#           axis.text.x = element_text(size = 20),
#           axis.title.x = element_text(size = 25),
#           axis.text.y = element_text(size = 20),
#           axis.title.y = element_text(size = 25),
#           legend.position = "none") +
#     scale_color_manual("Subtype",values=c(LumA = "#2176d5", LumB = "#34c6eb", Basal ="#c41b0e")) +
#     scale_x_continuous(limits = c(0.2, 0.8),breaks=seq(0, 0.8, 0.1)) +
#     scale_y_continuous(limits = c(0,0.6),breaks=seq(0, 0.6, 0.1))))


