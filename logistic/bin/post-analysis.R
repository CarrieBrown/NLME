# NLME Comparision Study
# Logistic Model - Simuation Metadata Analysis
# Carrie Brown - 3/21/20

rm(list=ls())

df <- 6
library(plyr)
library(tidyverse)
library(latex2exp)

source("../bin/nlme-library.R")

iml_parms <- read.csv("iml_varest_total.csv")
iml_parms$method <- "iml"
iml_yhat <- read.csv("iml_varest2_total.csv")
iml_yhat$method <- "iml"

nlm_temp <- read.csv("nlm_varest_total.csv")
nlm_temp <- nlm_temp %>% filter(!(substr(Parameter,1,3) == "log"))
nlm_parms <- data.frame(a = nlm_temp$Estimate[nlm_temp$Parameter == "alpha"])
nlm_parms$b <- nlm_temp$Estimate[nlm_temp$Parameter == "beta"]
nlm_parms$c <- nlm_temp$Estimate[nlm_temp$Parameter == "gamma"]
nlm_parms$d <- nlm_temp$Estimate[nlm_temp$Parameter == "delta"]
nlm_parms$se_a <- nlm_temp$StandardError[nlm_temp$Parameter == "alpha"]
nlm_parms$se_b <- nlm_temp$StandardError[nlm_temp$Parameter == "beta"]
nlm_parms$se_c <- nlm_temp$StandardError[nlm_temp$Parameter == "gamma"]
nlm_parms$se_d <- nlm_temp$StandardError[nlm_temp$Parameter == "delta"]

nlm_temp <- read.csv("nlm_varest2_total.csv")
nlm_yhat <- nlm_temp %>% filter(substr(Label, 1, 5) == "y-hat")
nlm_yhat$xi <- as.numeric(strsplit(as.character(nlm_yhat$Label), "=") %>% map_chr(2))
nlm_yhat <- nlm_yhat %>% select("xi", "Estimate", "StandardError")
colnames(nlm_yhat) <- c("xi", "yhat_xi", "se_yhat_xi")
nlm_yhat$method <- "nlm"

setwd("./meta_analysis_output")

var_temp <- data.frame(var_ai = nlm_temp %>% filter(substr(Label, 1, 2) == "ai") %>% select(Estimate),
                       var_bi = nlm_temp %>% filter(substr(Label, 1, 2) == "bi") %>% select(Estimate),
                       var_ci = nlm_temp %>% filter(substr(Label, 1, 2) == "ci") %>% select(Estimate),
                       var_di = nlm_temp %>% filter(substr(Label, 1, 2) == "di") %>% select(Estimate),
                       var_res = nlm_temp %>% filter(substr(Label, 1, 2) == "eu") %>% select(Estimate))
colnames(var_temp) <- c("var_ai", "var_bi", "var_ci", "var_di", "var_res")
nlm_parms <- cbind(nlm_parms, var_temp)
nlm_parms$method <- "nlm"

parms <- rbind(nlm_parms[,order(colnames(nlm_parms))], iml_parms[,order(colnames(iml_parms))])
yhat <- rbind(nlm_yhat, iml_yhat)
rm(nlm_temp, var_temp, nlm_parms, iml_parms, nlm_yhat, iml_yhat)

iml_converged <- nrow(parms %>% filter(method == "iml"))
nlm_converged <- nrow(parms %>% filter(method == "nlm"))

converged <- c(iml_converged, nlm_converged)
names(converged) <- c("iml_converged", "nlm_converged")
write.table(converged, "converged.csv", quote=FALSE)

# summary and plots of fixed effects

parm_list <- c("a", "b", "c", "d")
parm_label <- c("alpha", "beta", "gamma", "delta")
parm_file <- c("ahat", "bhat", "chat", "dhat")
parm_values <- c(10, 30, -7, 0.75)
n <- length(parm_list)
fixed_parm_capture <- data.frame(parameter = parm_list, 
                                 iml_captured = as.numeric(rep(NA, n)),  
                                 iml_high = as.numeric(rep(NA, n)),  
                                 iml_low = as.numeric(rep(NA, n)), 
                                 nlm_captured = as.numeric(rep(NA, n)),  
                                 nlm_high = as.numeric(rep(NA, n)),  
                                 nlm_low = as.numeric(rep(NA, n)))

for (i in 1:n){
  subset <- parms %>% select(ends_with(parm_list[i]), method)
  colnames(subset)[1] <- "x"
  newcol <- subset %>% select(starts_with("se_"))
  names(newcol) <- "se_x"
  subset <- cbind(subset, newcol)
  subset$method <- as.factor(subset$method)
  crit <- qt(0.975, df)
  subset <- subset %>% mutate(upper = x + crit*se_x, lower = x - crit*se_x)
  subset$capture <- as.factor("Y")
  levels(subset$capture) <- c("Y", "H", "L")
  subset$capture[which(subset$lower > parm_values[i])] <- "H"
  subset$capture[which(subset$upper < parm_values[i])] <- "L"
  
  fixed_parm_capture$iml_captured[i] <- nrow(subset %>% filter(method == "iml") %>% filter(capture == "Y"))
  fixed_parm_capture$iml_high[i] <- nrow(subset %>% filter(method == "iml") %>% filter(capture == "H"))
  fixed_parm_capture$iml_low[i] <- nrow(subset %>% filter(method == "iml") %>% filter(capture == "L"))
  fixed_parm_capture$nlm_captured[i] <- nrow(subset %>% filter(method == "nlm") %>% filter(capture == "Y"))
  fixed_parm_capture$nlm_high[i] <- nrow(subset %>% filter(method == "nlm") %>% filter(capture == "H"))
  fixed_parm_capture$nlm_low[i] <- nrow(subset %>% filter(method == "nlm") %>% filter(capture == "L"))
  
  title <- TeX(paste0("Distribution of $\\hat{\\",parm_label[i],"_i}$"))
  xlab <- TeX(paste0("$\\hat{\\",parm_label[i],"_i$"))
  plot_parm_dist(subset, parm_list[i], title, xlab, parm_file[i], parm_values[i])

  plot_parm_ci(subset, parm_list[i], parm_file[i], parm_values[i])
  
  temp <- as.data.frame(summary(subset %>% filter(method == "iml") %>% select(x)))
  newrow <- list(c(parm_list[i], "iml", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
  names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
  newrow <- unlist(newrow)
  
  if (i == 1)
    fixed_parm_summary <- data.frame(t(newrow), stringsAsFactors = FALSE)
  else
    fixed_parm_summary <- rbind(fixed_parm_summary, newrow)
  
  temp <- as.data.frame(summary(subset %>% filter(method == "nlm") %>% select(x)))
  newrow <- list(c(parm_list[i], "nlm", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
  names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
  newrow <- unlist(newrow)
  
  fixed_parm_summary <- rbind(fixed_parm_summary, newrow)

}

fixed_parm_summary
write.csv(fixed_parm_summary, "parameter_analysis_summary.csv",quote=FALSE, row.names = FALSE)

fixed_parm_capture
write.csv(fixed_parm_capture, "parameter_analysis_capture.csv",quote=FALSE, row.names = FALSE)

# plot distributions of variance estimates

# random effects
parm_list <- c("var_ai", "var_bi", "var_ci", "var_di")
parm_label <- c("alpha", "beta", "gamma", "delta")
parm_file <- parm_list
parm_values <- c(0.5, 2.5, 0.05, 0.005)

n <- length(parm_list)
random_var_summary <- data.frame()

for (i in 1:n){
  subset <- parms %>% select(ends_with(parm_list[i]), method)
  colnames(subset)[1] <- "x"
  subset$method <- as.factor(subset$method)
  
  title <- TeX(paste0("Distribution of $\\widehat{\\Var(\\",parm_label[i],")}$"))
  xlab <- TeX(paste0("$\\widehat{Var(\\",parm_label[i],")$"))
  plot_parm_dist(subset, parm_list[i], title, xlab, parm_file[i], parm_values[i])
  
  temp <- as.data.frame(summary(subset %>% filter(method == "iml") %>% select(x)))
  newrow <- list(c(parm_list[i], "iml", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
  names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
  newrow <- unlist(newrow)
  
  if (i == 1)
    random_var_summary <- data.frame(t(newrow), stringsAsFactors = FALSE)
  else
    random_var_summary <- rbind(random_var_summary, newrow)
    
  temp <- as.data.frame(summary(subset %>% filter(method == "nlm") %>% select(x)))
  newrow <- list(c(parm_list[i], "nlm", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
  names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
  newrow <- unlist(newrow)
  
  random_var_summary <- rbind(random_var_summary, newrow)
}

# residual

subset <- parms %>% select(var_res, method)
colnames(subset)[1] <- "x"
subset$method <- as.factor(subset$method)

title <- TeX(paste0("Distribution of $\\widehat{\\Var(eu_i)}$"))
xlab <- TeX(paste0("$\\widehat{Var(\\eu_i)$"))
plot_parm_dist(subset, "var_eu", title, xlab, "var_eu", 0.5)

temp <- as.data.frame(summary(subset %>% filter(method == "iml") %>% select(x)))
newrow <- list(c("var_eu", "iml", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
newrow <- unlist(newrow)
random_var_summary <- rbind(random_var_summary, newrow)

temp <- as.data.frame(summary(subset %>% filter(method == "nlm") %>% select(x)))
newrow <- list(c("var_eu", "nlm", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
newrow <- unlist(newrow)
random_var_summary <- rbind(random_var_summary, newrow)

random_var_summary
write.csv(random_var_summary, "variance_analysis_summary.csv",quote=FALSE, row.names = FALSE)

# Yhat analysis

parm_list <- unique(yhat$xi)
parm_label <- as.character(unique(yhat$xi))
parm_file <- paste0("yhat_",parm_label)
parm_values <- 10 + (30 / (1 + exp(7 - (0.75 * parm_list))))
n <- length(parm_list)
yhat_capture <- data.frame(parameter = parm_list, 
                                 iml_captured = as.numeric(rep(NA, n)),  
                                 iml_high = as.numeric(rep(NA, n)),  
                                 iml_low = as.numeric(rep(NA, n)), 
                                 nlm_captured = as.numeric(rep(NA, n)),  
                                 nlm_high = as.numeric(rep(NA, n)),  
                                 nlm_low = as.numeric(rep(NA, n)))

for (i in 1:n){
  subset <- yhat %>% filter(xi == parm_list[i]) %>% select(yhat_xi, se_yhat_xi, method)
  colnames(subset) <- c("x", "se_x", "method")
  subset$method <- as.factor(subset$method)
  subset <- subset %>% mutate(upper = x + crit*se_x, lower = x - crit*se_x)
  subset$capture <- as.factor("Y")
  levels(subset$capture) <- c("Y", "H", "L")
  subset$capture[which(subset$lower > parm_values[i])] <- "H"
  subset$capture[which(subset$upper < parm_values[i])] <- "L"
  
  yhat_capture$iml_captured[i] <- nrow(subset %>% filter(method == "iml") %>% filter(capture == "Y"))
  yhat_capture$iml_high[i] <- nrow(subset %>% filter(method == "iml") %>% filter(capture == "H"))
  yhat_capture$iml_low[i] <- nrow(subset %>% filter(method == "iml") %>% filter(capture == "L"))
  yhat_capture$nlm_captured[i] <- nrow(subset %>% filter(method == "nlm") %>% filter(capture == "Y"))
  yhat_capture$nlm_high[i] <- nrow(subset %>% filter(method == "nlm") %>% filter(capture == "H"))
  yhat_capture$nlm_low[i] <- nrow(subset %>% filter(method == "nlm") %>% filter(capture == "L"))
  
  title <- TeX(paste0("Distribution of $\\widehat{\\Y_{x=",parm_label[i],"}}$"))
  xlab <- TeX(paste0("$\\widehat{\\Y_{x=",parm_label[i],"}}$"))
  plot_parm_dist(subset, parm_list[i], title, xlab, parm_file[i], parm_values[i])
  
  plot_parm_ci(subset, parm_list[i], parm_file[i], parm_values[i])
  
  temp <- as.data.frame(summary(subset %>% filter(method == "iml") %>% select(x)))
  newrow <- list(c(parm_list[i], "iml", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
  names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
  newrow <- unlist(newrow)
  
  if (i == 1)
    yhat_summary <- data.frame(t(newrow), stringsAsFactors = FALSE)
  else
    yhat_summary <- rbind(yhat_summary, newrow)
  
  temp <- as.data.frame(summary(subset %>% filter(method == "nlm") %>% select(x)))
  newrow <- list(c(parm_list[i], "nlm", as.numeric(strsplit(as.character(temp$Freq), ":") %>% map_chr(2))))
  names(newrow[[1]]) <- c("parm", "method", "min", "firstq", "median", "mean", "thirdq", "max")
  newrow <- unlist(newrow)
   
  yhat_summary <- rbind(yhat_summary, newrow)
}

yhat_summary
write.csv(yhat_summary, "yhat_analysis_summary.csv",quote=FALSE, row.names = FALSE)

yhat_capture
write.csv(yhat_capture, "yhat_analysis_capture.csv",quote=FALSE, row.names = FALSE)


