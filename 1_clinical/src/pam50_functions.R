# Functions for PAM50 project

# loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
}

# define "not in" operator 
'%!in%' <- function(x,y)!('%in%'(x,y))

# simple KM plot function
# function that creates KM plot for specified OM (2 or 3 lines)
KM.plot <- function(data,fit,OM,colors,title,group.variable,legend.labels,legend.title = "majoritySecondBestClass") { #
  
  legend.labs = if(length(legend.labels)==2) {
    c(paste(legend.labels[1]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[1],")",sep = ""),
      paste(legend.labels[2]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[2],")",sep = ""))
  } else if(length(legend.labels)==3) {
    c(paste(legend.labels[1]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[1],")",sep = ""),
   paste(legend.labels[2]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[2],")",sep = ""),
   paste(legend.labels[3]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[3],")",sep = ""))
  } else if(length(legend.labels)==4) {
    c(paste(legend.labels[1]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[1],")",sep = ""),
      paste(legend.labels[2]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[2],")",sep = ""),
      paste(legend.labels[3]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[3],")",sep = ""),
      paste(legend.labels[4]," (",table(data[!is.na(data[[OM]]),][[group.variable]])[4],")",sep = ""))
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
    xlim = c(0,max(data[[OM]][is.finite(data[[OM]])])),         
    xlab = paste(OM," (years)", sep = ""),
    ylab = paste(OM," event probability", sep = ""), 
    ylim = c(0,1),
    color = "strata",
    palette = colors, 
    legend = c(0.8,0.2),
    ggtheme = theme(legend.title = element_text(size=25), #20
                    legend.key.size = unit(0.5,"cm"), 
                    legend.text = element_text(size = 25), #20
                    axis.text.x = element_text(size = 25), #20
                    axis.title.x = element_text(size = 30), #25
                    axis.text.y = element_text(size = 25), #20
                    axis.title.y = element_text(size = 30),
                    plot.title = element_text(size=30)),
    title = title,
    legend.title = legend.title,
    # HERE HOW TO HANDLE DIFFERENT NUMBERS OF LABELS
    legend.labs = legend.labs,
    break.x.by = 1, # break X axis in time intervals of x (what is nicest here? maybe 365)
    break.y.by = 0.1)
  
  return(plot)
  
}

# function that created univariate Cox proportional hazards model forest plot
unicox <- function(data,surv,title) {
  
  # Model construction
  main.pam50 <- coxph(surv~PAM50, data=data)
  
  # Result
  print(summary(main.pam50))
  
  # forest 
  ggforest(main.pam50,fontsize = 2,data=data,main=title)
}

# function that created multivariate Cox proportional hazards model forest plot
mvcox <- function(data,surv,title) {
  
  # Model construction 
  # parameters to incl: PAM50, Age, Grade, TumSize
  main.all <- coxph(surv~PAM50+Age+LN+TumSize+Grade, data=data) 
  
  # Result
  print(summary(main.all))
  
  # Plot forest 
  ggforest(main.all,fontsize = 2,cpositions = c(0.01,0.13,0.35),data=data,main=title)
}