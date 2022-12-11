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

# plot list
plot.list <- list()

#######################################################################
# Define groups
#######################################################################

pam50.subtypes <- unique(loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData")$NCN_PAM50)

clinical.groups <- setNames(as.data.frame(rbind(
  All = c(NA,NA), 
  TNBC = c("Negative","Negative"), 
  ERpHER2p = c("Positive","Positive"), 
  ERnHER2p = c("Negative","Positive"), 
  ERpHER2n = c("Positive","Negative")), 
  stringsAsFactors = FALSE), 
  c("ER","HER2"))

#######################################################################
# PAM50 centroid correlation - NOT CHECK!!!!!!
#######################################################################

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData")

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
    pam50.data.mod <- melt(pam50.data, measure.vars=colnames(pam50.data[grepl("mean",colnames(pam50.data))])) # only work if other columns dont caintain "mean" substring
    count.df <- count(pam50.data.mod,variable) %>% mutate(y_pos = min(pam50.data.mod$value))
    # plot
    plot.list <- append(plot.list,list(
      ggplot(pam50.data.mod) + 
        geom_boxplot(aes(x=variable, y=value,fill=as.factor(variable)),alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("PAM50 Centroid") +
        ylab("Correlation") +
        ggtitle(paste(
          "PAM50 centroid correlation - clin: ",clin.group,"; pam50: ",pam50.type," ","(n=",unname(table(pam50.data$NCN_PAM50)[1]),")",sep="")) +
        geom_text(data = count.df, aes(x=variable, y = y_pos, label = paste("n=",n,sep="")),vjust=1.8,size=8) + 
        scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
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
# Crosstable with 2nd best PAM50 match for all samples - CHECK
#######################################################################

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData")

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
# best correlation matches - CHECK
#######################################################################

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData")

# filter out samples where best=2nd best and calc diff.
anno <- anno[which(
  anno$majorityNearestClass!=anno$majoritySecondBestClass),] %>% 
  rowwise() %>% 
  mutate(Diff = abs(
    abs(get(paste("mean",majorityNearestClass,sep=""))) - abs(get(paste("mean",majoritySecondBestClass,sep="")))
    ))

for (i in 1:nrow(clinical.groups)) { 
  # select subgroup data
  clin.group <- row.names(clinical.groups)[i] 

  if (clin.group == "All") {
    sg.data <- anno
  } else {
    sg.data <- anno[which(
      anno$HER2 == clinical.groups[clin.group ,"HER2"] & 
        anno$ER == clinical.groups[clin.group ,"ER"]),]
  }
  count.df <- count(sg.data,majorityNearestClass) %>% 
    mutate(y_pos = min(sg.data$Diff))
  
  # plot
  plot.list <- append(plot.list,list(
    ggplot(sg.data) + 
    geom_boxplot(aes(x=majorityNearestClass, y=Diff,fill=as.factor(majorityNearestClass)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 Class") +
    ylab("Correlation diff (1st-2nd)") +
    ggtitle(paste("PAM50 class distinctiveness - clin: ",clin.group,sep="")) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
      geom_text(data = count.df, aes(x=majorityNearestClass, y = y_pos, label = paste("n=",n,sep="")),vjust=1.8,size=8) +
    scale_fill_manual(values= c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb", Basal ="#c41b0e", Normal="#64c25d"))))
}

#######################################################################
# PAM50 Class specific distinctiveness: difference between the best and second best correlation matches - CHECK
#######################################################################

# plot list 
csd.list <- list()

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData")

# filter out samples where best=2nd best and calc diff.
anno <- anno[which(
  anno$majorityNearestClass!=anno$majoritySecondBestClass),] %>% 
  rowwise() %>% 
  mutate(Diff = abs(
    abs(get(paste("mean",majorityNearestClass,sep=""))) - abs(get(paste("mean",majoritySecondBestClass,sep="")))
    ))

for (i in 1:nrow(clinical.groups)) { 
  # select subgroup data
  clin.group <- row.names(clinical.groups)[i] 
  
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
    count.df <- count(pam50.data,majoritySecondBestClass) %>% 
      mutate(y_pos = min(pam50.data$Diff))
    
    # plot
    plot.list <- append(plot.list,list(
      ggplot(pam50.data) + 
        geom_boxplot(aes(x=majoritySecondBestClass, y=Diff,fill=as.factor(majoritySecondBestClass)),alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("majoritySecondBestClass") +
        ylab("Correlation diff (1st-2nd)") +
        ggtitle(paste("PAM50 class-specific distinctiveness - clin: ",clin.group,"; pam50: ",pam50.type,sep="")) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        theme(plot.title = element_text(size = 25),
              axis.text.x = element_text(size = 20),
              axis.title.x = element_text(size = 25),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 25),
              legend.position = "none") +
        geom_text(data = count.df, aes(x=majoritySecondBestClass, y = y_pos, label = paste("n=",n,sep="")),vjust=1.8,size=8) + #
        scale_fill_manual(values= c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb", Basal ="#c41b0e", Normal="#64c25d"))))
  }
}

#######################################################################
# SA for second best pam50 class - CHECK  
#######################################################################

# plot list 
sa.list <- list()

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

#######################################################################
# SA for distinctiveness in LumA and LumB
#######################################################################

# load data
anno <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData") %>% 
  filter(EvalGroup_ERpHER2nLNn50 == "ERpHER2nLNn_Endo_50") %>% 
  filter(majorityNearestClass != majoritySecondBestClass) %>% 
  rowwise() %>% 
  mutate(Diff = abs(
    abs(get(paste("mean",majorityNearestClass,sep=""))) - abs(get(paste("mean",majoritySecondBestClass,sep="")))
  ))

# LumA
# KM plot: LumA-Normal (distinct); LumA-Normal (non-distinct); LumA-LumB (distinct); LumA-LumB (non-distinct)

median.LumA.LumB <- median(
  anno[which(anno$majorityNearestClass =="LumA" &
               anno$majoritySecondBestClass=="LumB"),]$Diff)
median.LumA.Normal <- median(
  anno[which(anno$majorityNearestClass =="LumA" &
               anno$majoritySecondBestClass=="Normal"),]$Diff)
luma.data <- anno %>% 
  filter(majorityNearestClass == "LumA") %>% 
  filter(majoritySecondBestClass %in% c("LumB","Normal")) %>% 
  mutate(Group = case_when(
    majoritySecondBestClass=="LumB" & Diff >= median.LumA.LumB ~ "LumB-Distinct",
    majoritySecondBestClass=="LumB" & Diff <= median.LumA.LumB ~ "LumB-NonDistinct",
    majoritySecondBestClass=="Normal" & Diff >= median.LumA.Normal ~ "Normal-Distinct",
    majoritySecondBestClass=="Normal" & Diff <= median.LumA.Normal ~ "Normal-NonDistinct"
                           ))

# plot KM
for (k in 1:length(outcome.measures)) {
  OM <- outcome.measures[k]
  OM.bin <- paste(OM, "bin",sep="")    
  # surv object
  data.surv <- Surv(luma.data[[OM]], luma.data[[OM.bin]]) 
  # fit
  fit <- survminer::surv_fit(data.surv~Group, data=luma.data, conf.type="log-log")
  # plot
  plot.list <- append(plot.list,list(
    KM.plot(luma.data,fit,OM,
            colors = c("#2b8cbe","#7bccc4","#238443","#78c679"),
            title = paste(OM, " in Luminal A 2nd-best/distinctiveness subgroups (ERpHER2nLNn_Endo_50_LumA)", sep=""),
            group.variable = "Group",
            legend.labels = names(table(luma.data$Group)),
            legend.title = "Subgroups")
    )) 
}

# Boxplots: Age; Tum Size (with these groups)

# plot
# age
count.df <- count(luma.data,Group) %>% 
  mutate(y_pos = min(luma.data$Age))
plot.list <- append(plot.list,list(
  ggplot(luma.data) + 
    geom_boxplot(aes(x=Group, y=Age,fill=as.factor(Group)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("Second-best PAM50 class - LumA Distinctiveness") +
    ylab("Age (years)") +
    ggtitle("Age in ERpHER2nLNn_Endo_50; NCN_PAM50: Luminal A") +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    geom_text(data = count.df, aes(x=Group, y = y_pos, label = paste("n=",n,sep="")),vjust=1.6,size=8) + #
    scale_fill_manual(values= c("#2b8cbe","#7bccc4","#238443","#78c679"))))

# size Size_mm
count.df <- count(luma.data,Group) %>% 
  mutate(y_pos = min(luma.data$Size_mm))
plot.list <- append(plot.list,list(
  ggplot(luma.data) + 
    geom_boxplot(aes(x=Group, y=Size_mm,fill=as.factor(Group)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("Second-best PAM50 class - LumA Distinctiveness") +
    ylab("Size (mm)") +
    ggtitle("Tumor size in ERpHER2nLNn_Endo_50; NCN_PAM50: Luminal A") +
    scale_y_continuous(breaks = seq(0, 70, 10),limits=c(0,70)) +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    geom_text(data = count.df, aes(x=Group, y = y_pos, label = paste("n=",n,sep="")),vjust=1.6,size=8) + #
    scale_fill_manual(values= c("#2b8cbe","#7bccc4","#238443","#78c679"))))


# LumB
# KM plot: LumB-Her2 (distinct); LumB-Her2 (non-distinct); LumB-LumA (distinct); LumB-LumA (non-distinct)

median.LumB.LumA <- median(
  anno[which(anno$majorityNearestClass =="LumB" &
               anno$majoritySecondBestClass=="LumA"),]$Diff)
median.LumB.Her2 <- median(
  anno[which(anno$majorityNearestClass =="LumB" &
               anno$majoritySecondBestClass=="Her2"),]$Diff)
lumb.data <- anno %>% 
  filter(majorityNearestClass == "LumB") %>% 
  filter(majoritySecondBestClass %in% c("LumA","Her2")) %>% 
  mutate(Group = case_when(
    majoritySecondBestClass=="LumA" & Diff >= median.LumB.LumA ~ "LumA-Distinct",
    majoritySecondBestClass=="LumA" & Diff <= median.LumB.LumA ~ "LumA-NonDistinct",
    majoritySecondBestClass=="Her2" & Diff >= median.LumB.Her2 ~ "Her2-Distinct",
    majoritySecondBestClass=="Her2" & Diff <= median.LumB.Her2 ~ "Her2-NonDistinct"
  ))

# plot KM
for (k in 1:length(outcome.measures)) {
  OM <- outcome.measures[k]
  OM.bin <- paste(OM, "bin",sep="")    
  # surv object
  data.surv <- Surv(lumb.data[[OM]], lumb.data[[OM.bin]]) 
  # fit
  fit <- survminer::surv_fit(data.surv~Group, data=lumb.data, conf.type="log-log")
  # plot
  plot.list <- append(plot.list,list(
    KM.plot(lumb.data,fit,OM,
            colors = c("#ae017e","#f768a1","#225ea8","#41b6c4"), 
            title = paste(OM, " in Luminal B 2nd-best/distinctiveness subgroups (ERpHER2nLNn_Endo_50_LumB)", sep=""),
            group.variable = "Group",
            legend.labels = names(table(lumb.data$Group)),
            legend.title = "Subgroups")
    )) 
}

# Boxplots: Age; Tum Size (with these groups)

# plot
# age
count.df <- count(lumb.data,Group) %>% 
  mutate(y_pos = min(lumb.data$Age))
plot.list <- append(plot.list,list(
  ggplot(lumb.data) + 
    geom_boxplot(aes(x=Group, y=Age,fill=as.factor(Group)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("Second-best PAM50 class - LumB Distinctiveness") +
    ylab("Age (years)") +
    ggtitle("Age in ERpHER2nLNn_Endo_50; NCN_PAM50: Luminal B") +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    geom_text(data = count.df, aes(x=Group, y = y_pos, label = paste("n=",n,sep="")),vjust=1.6,size=8) + #
    scale_fill_manual(values= c("#ae017e","#f768a1","#225ea8","#41b6c4"))))

# size Size_mm
count.df <- count(lumb.data,Group) %>% 
  mutate(y_pos = min(lumb.data$Size_mm))
plot.list <- append(plot.list,list(
  ggplot(lumb.data) + 
    geom_boxplot(aes(x=Group, y=Size_mm,fill=as.factor(Group)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("Second-best PAM50 class - LumB Distinctiveness") +
    ylab("Size (mm)") +
    ggtitle("Tumor size in ERpHER2nLNn_Endo_50; NCN_PAM50: Luminal B") +
    scale_y_continuous(breaks = seq(0, 70, 10),limits=c(0,70)) +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    geom_text(data = count.df, aes(x=Group, y = y_pos, label = paste("n=",n,sep="")),vjust=1.6,size=8) + #
    scale_fill_manual(values= c("#ae017e","#f768a1","#225ea8","#41b6c4"))))

# print to pdf
pdf(file = paste(output.path,"PAM50_analyses.pdf", sep=""),
    onefile = TRUE, width = 21.0, height = 14.8)

for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}


dev.off()
