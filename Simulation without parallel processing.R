library(mice)
library(reshape2)
library(ggplot2)
library(stringr)
theme_set(theme_classic())
set.seed(14850)

DAG_xy_collider <- function(N, confounding = F) {
  V1 <- rnorm(N)
  u.V2 <- rnorm(N)
  V2 <- u.V2 + rnorm(N)
  u.c <- rnorm(N)
  collider <- V1 + V2 + u.c + rnorm(N)
  
  df <- as.data.frame(cbind(V2,V1,collider,u.V2,u.c))
  return(df)
  if(confounding) {print(1)}
}

DAG_xz_collider <- function(N, confounding = F) {
  u.z <- rnorm(N)
  z <- u.z + rnorm(N)
  V1 <- rnorm(N)
  if (confounding) { V1 <- z + V1 }
  V2 <- z + rnorm(N)
  u.c <- rnorm(N)
  collider <- u.c + V1 + z + rnorm(N)
  
  df <- as.data.frame(cbind(V1,V2,z,collider,u.z,u.c))
  return(df)
}

DAG_M_bias <- function(N, confounding = F) {
  u.z <- rnorm(N)
  z <- u.z + rnorm(N)
  u.z_2 <- rnorm(N)
  z_2 <- u.z_2 + rnorm(N)
  u.c <- rnorm(N)
  collider <- u.c + z + z_2 + rnorm(N)
  V1 <- z + rnorm(N)
  V2 <- z_2 + rnorm(N)
  
  df <- as.data.frame(cbind(V1,V2,z,z_2,collider,u.z,u.z_2,u.c))
  return(df)
  if(confounding) {print(1)}
}

summary_lm <- function(formula, data) {
  return(summary(lm(as.formula(formula), data=data))$coefficients)
}
summary_lm_mice <- function(formula, data) {
  return(summary(pool(with(data, lm(as.formula(formula))))))
}

simulation_missingness <- function(n_replicates, N, DAG, formula = "V2 ~ V1 + z", missing_var = "z", 
                                     missingness_mechanism = c("V1", "u.c"), confounding = F){
  # set up the data 
  dat.full <- DAG(N, confounding)
  
  # create missingness in [missing_var]
  dat <- dat.full
  miss.indicator <- rowSums(dat[missingness_mechanism]) # take the sum of the variables that cause missingness
  miss.indicator <- pnorm(miss.indicator) # get a probability of missingness from the magnitude of said sum
  miss.indicator <- rbinom(N,1,miss.indicator) # binomial for missingness or not from said probability
  dat[miss.indicator==1,missing_var] <- NA
  
  #c Create different datasets before MICE leaving out certain columns
  dat.no_uz <- dat[!names(dat) %in% "u.z"]
  dat.no_uc <- dat[!names(dat) %in% "u.c"]
  dat.no_uc_or_uz <- dat[!names(dat) %in% c("u.c", "u.z")]
  dat.nocol <- dat[!names(dat) %in% "collider"]
  dat.nocol_nouz <- dat[!names(dat) %in% c("collider","u.z")]
  dat.nocol_nouc <- dat[!names(dat) %in% c("collider", "u.c")]
  dat.nocol_nouc_nouz <- dat[!names(dat) %in% c("collider", "u.c", "u.z")]
  
  # impute the data
  mi <- mice(dat, printFlag = FALSE)         
  mi.no_uz <- mice(dat.no_uz, printFlag = FALSE) 
  mi.no_uc <- mice(dat.no_uc, printFlag = FALSE)
  mi.no_uc_or_uz <- mice(dat.no_uc_or_uz, printFlag = FALSE)
  mi.nocol <- mice(dat.nocol, printFlag = FALSE)
  mi.nocol_nouz <- mice(dat.nocol_nouz, printFlag = FALSE)
  mi.nocol_nouc <- mice(dat.nocol_nouc, printFlag = FALSE)
  mi.nocol_nouc_nouz <- mice(dat.nocol_nouc_nouz, printFlag = FALSE)
  
  # Perform regressions
  res.full <- summary_lm(formula, dat.full) # if all data were observed
  res.ld <- summary_lm(formula, dat) # with listwise deletion
  
  res.mi <- summary_lm_mice(formula, mi) # with all nodes
  res.mi.no_uz <- summary_lm_mice(formula, mi.no_uz) # without u.z
  res.mi.no_uc <- summary_lm_mice(formula, mi.no_uc) # without u.c
  res.mi.no_uc_or_uz <- summary_lm_mice(formula, mi.no_uc_or_uz) # without u.z or u.c
  
  if (str_detect(formula, "collider")) {formula <- str_remove(formula, "\\+collider")}
  res.mi.nocol <- summary_lm_mice(formula, mi.nocol) # without collider
  res.mi.nocol_nouz <- summary_lm_mice(formula, mi.nocol_nouz) # without collider or u.z
  res.mi.nocol_nouc <- summary_lm_mice(formula, mi.nocol_nouc) # without collider or u.c
  res.mi.nocol_nouc_nouz <- summary_lm_mice(formula, mi.nocol_nouc_nouz) # without collider, u.z or u.c
  
  # Get results
  r <- 2
  ci <- c(1,2,4,6)
  outputs <- cbind(res.full[r,ci[1]],
                   res.ld[r,ci[1]],
                   res.mi[r,ci[2]],
                   res.mi.no_uz[r,ci[2]],
                   res.mi.no_uc[r,ci[2]],
                   res.mi.no_uc_or_uz[r,ci[2]],
                   res.mi.nocol[r,ci[2]],
                   res.mi.nocol_nouz[r,ci[2]],
                   res.mi.nocol_nouc[r,ci[2]],
                   res.mi.nocol_nouc_nouz[r,ci[2]],
                   res.full[r,ci[3]],
                   res.ld[r,ci[3]],
                   res.mi[r,ci[4]],
                   res.mi.no_uz[r,ci[4]],
                   res.mi.no_uc[r,ci[4]],
                   res.mi.no_uc_or_uz[r,ci[4]],
                   res.mi.nocol[r,ci[4]],
                   res.mi.nocol_nouz[r,ci[4]],
                   res.mi.nocol_nouc[r,ci[4]],
                   res.mi.nocol_nouc_nouz[r,ci[4]])
  cat(i, "out of", n_replicates," \r")
  i <<- i+1
  return(outputs)
}

# Function to replicate the simulation n number of times
repeats <- function(n_replicates, N, input_DAG, formula, missing_var, missingness_mechanism, confounding = F) {
  t1 <- paste0("DAG: ", as.character(substitute(input_DAG)))
  t2 <- paste0("Formula: ", formula)
  t3 <- paste0("Missingness in: ", missing_var)
  t4 <- paste0("Causes of missingness: ", paste(missingness_mechanism, collapse=", "))
  print(c(t1,t2,t3,t4))
  
  i <<- 1
  
  input_DAG <- input_DAG # necessary step as using a function from argument directly into an argument causes error
  sims <- data.frame(t(matrix(replicate(n_replicates, simulation_missingness(n_replicates, N, input_DAG, formula, missing_var, missingness_mechanism, confounding)),nrow=20)))
  
  
  sims_estimate <- sims[1:10]
  sims_pval <- sims[11:20]
  
  sims_estimate <- cbind(1:nrow(sims_estimate), sims_estimate)
  sims_pval <- cbind(1:nrow(sims_pval), sims_pval)
  
  names(sims_estimate) <- names(sims_pval) <- c("Iter", "Full","LD","MI","MI-no u.z", 
                                                "MI-no u.c","MI-no  u.c or u.z", 
                                                "MI-no c", "MI-no c or u.z", 
                                                "MI-no c or u.c", "MI-no c or u.c or u.z")
  
  # If DAG_xy_collider was used, remove all estimators which excluded z or u.z 
  # (as these variables don't exist anyway)
  if (identical(input_DAG, DAG_xy_collider)) {
    col_rm <- grepl("z", colnames(sims_estimate))
    sims_estimate <- sims_estimate[!col_rm]
    sims_pval <- sims_pval[!col_rm]
  }
  
  return(list(Estimates = sims_estimate, P_vals = sims_pval, Description = c(t1,t2,t3,t4)))
}

# Function to create plots from sim results
sim_plot <- function(result, start = 2, end = ncol(result$Estimates)) {
  info <- result$Description
  result <- result$Estimates[c(1,start:end)]
  to.plot.b <- melt(result, id.vars="Iter", variable.name="Estimator", value.name="Estimate")
  ggplot(to.plot.b, aes(x=Estimate, color=Estimator)) + geom_density() + ggtitle(paste(info, collapse=".\n"))
}

# Function to print means estimates and false positive rates
sim_print_values <- function(input) {
  ests <- input$Estimates
  pvals <- input$P_vals
  print("Estimates:")
  print(round(colMeans(ests[2:ncol(ests)]),2))
  print("False positive rates:")
  print(colMeans(pvals[2:ncol(pvals)]<.05))
}

####
t1.1 <- repeats(250, 1000, DAG_xy_collider, "V2~V1", "V2", c("V1", "u.c"))
t1.2 <- repeats(250, 1000, DAG_xy_collider, "V2~V1", "V2", c("V1", "collider"))
t1.3 <- repeats(250, 1000, DAG_xy_collider, "V2~V1", "V2", c("V1", "u.V2"))

t1.4 <- repeats(250, 1000, DAG_xy_collider, "V1~V2", "V1", c("u.V2", "u.c"))
t1.5 <- repeats(250, 1000, DAG_xy_collider, "V1~V2", "V1", c("u.V2", "collider"))
t1.6 <- repeats(250, 1000, DAG_xy_collider, "V1~V2", "V1", c("u.V2"))

t1.7 <- repeats(250, 1000, DAG_xy_collider, "V1~V2", "V2", c("u.V2", "u.c"))
t1.8 <- repeats(250, 1000, DAG_xy_collider, "V1~V2", "V2", c("u.V2", "collider"))
t1.9 <- repeats(250, 1000, DAG_xy_collider, "V1~V2", "V2", c("u.V2"))



####
t2.1 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "z", c("V1", "u.c"))
t2.2 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "z", c("V1", "collider"))
t2.3 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "z", c("V1", "u.z"))

t2.4 <- repeats(250, 1000, DAG_xz_collider, "V1~V2", "z", c("u.z", "u.c"))
t2.5 <- repeats(250, 1000, DAG_xz_collider, "V1~V2", "z", c("u.z", "collider"))
t2.6 <- repeats(250, 1000, DAG_xz_collider, "V1~V2", "z", c("u.z"))

t2.7 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "V2", c("V1", "u.c"))
t2.8 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "V2", c("V1", "collider"))
t2.9 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "V2", c("V1", "u.z"))

t2.11 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "z", c("V1", "u.c"), TRUE)
t2.12 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "z", c("V1", "collider"), TRUE)
t2.13 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "z", c("V1", "u.z"), TRUE)

t2.14 <- repeats(250, 1000, DAG_xz_collider, "V1~V2", "z", c("u.z", "u.c"), TRUE)
t2.15 <- repeats(250, 1000, DAG_xz_collider, "V1~V2", "z", c("u.z", "collider"), TRUE)
t2.16 <- repeats(250, 1000, DAG_xz_collider, "V1~V2", "z", c("u.z"), TRUE)

t2.17 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "V2", c("V1", "u.c"), TRUE)
t2.18 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "V2", c("V1", "collider"), TRUE)
t2.19 <- repeats(250, 1000, DAG_xz_collider, "V2~V1", "V2", c("V1", "u.z"), TRUE)


####
t3.1 <- repeats(250, 1000, DAG_M_bias, "V2~V1", "V2", c("V1", "u.c"))
t3.2 <- repeats(250, 1000, DAG_M_bias, "V2~V1", "V2", c("z"))
t3.3 <- repeats(250, 1000, DAG_M_bias, "V2~V1", "V2", c("u.z", "u.z_2"))

t3.4 <- repeats(250, 1000, DAG_M_bias, "V2~V1+collider+z_2", "z", c("V1", "u.c"))
t3.5 <- repeats(250, 1000, DAG_M_bias, "V2~V1+collider+z_2", "z", c("u.z","u.c"))
t3.6 <- repeats(250, 1000, DAG_M_bias, "V2~V1+collider+z_2", "z", c("V1", "u.z"))

t3.7 <- repeats(250, 1000, DAG_M_bias, "V1~V2+collider+z_2", "z", c("V2", "u.c"))
t3.8 <- repeats(250, 1000, DAG_M_bias, "V1~V2+collider+z_2", "z", c("u.z","u.c"))
t3.9s <- repeats(250, 1000, DAG_M_bias, "V1~V2+collider+z_2", "z", c("V2", "u.z"))


remove_columns <- function(input) {
  col_rm <- grepl("no c", colnames(input$Estimates))
  input$Estimates <- input$Estimates[!col_rm]
  input$P_vals <- input$P_vals[!col_rm]
  return(input)
}

t3.4_short <- remove_columns(t3.4)
t3.5_short <- remove_columns(t3.5)
t3.6_short <- remove_columns(t3.6)
t3.7_short <- remove_columns(t3.7)
t3.8_short <- remove_columns(t3.8)
t3.9_short <- remove_columns(t3.9)



####
sim_plot(t3.9_short)
sim_print_values(t3.9_short)