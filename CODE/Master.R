# Title     : Machine Learning Project Simulation Script
# Objective : Replicate simulations presented in Carrasco & Rossi (2016)
# Created by: Nikita Marini & Luca Poll
# Created on: 4/13/2021


## =============================================================================
## ==  Set environment      ====================================================
## =============================================================================

library(tidyverse)

# define paths
#setwd("..")
scripts <- paste0(getwd(), "/CODE/")
datapath <- paste0(getwd(), "/DATA/")
simoutpath <- paste0(getwd(), "/OUTPUT/SIMULATION/")

# load functions
source(paste0(scripts, "Simulation.R"))  # for DGP 1-6 Simulation function
source(paste0(scripts, "Estimators.R"))  # functions to get alpha range and estimators




## =============================================================================
## ==  Monte Carlo Simulation Function   =======================================
## =============================================================================

monte_carlo <- function(sample, estimator, size){

  # setup
  output <- data.frame(Estimator = character(), DGP = numeric(), IC = character(),
                       alpha = numeric(),
                       tuningparam = character(), tp = numeric(),
                       DoF = numeric(), MSE = numeric(),
                       RMSFE_rec = numeric(), RMSFE_roll = numeric(),
                       stringsAsFactors=FALSE)

  # assign feasible information criteria
  if (estimator %in% c("ridge", "LF")){
    ICs <- c("Mallows", "GCV")
  } else if (estimator == "PLS"){
    ICs <- "LOOCV"
  } else if (estimator == "PC"){
    ICs <- c("GCV", "Mallows", "AIC", "BIC", "AIC_sig", "BIC_sig", "Bai_Ng_pc", "Bai_Ng_ic")
  }

  #loop through simulation size (100 times)
  for (m in seq_len(size)){
    print(paste("m =", m))

    # 1) Generate data
    if (sample == "large"){
      data_full <- sampler()
    } else if (sample == "small"){
      data_full <- sampler(small_sample = T)
    }

    # loop through DGPs (6 times)
    for (d in seq_along(data_full)){
      #print(paste("d =", d))

      data <- data_full[[d]]

      # 2) get alpha range
      alpha <- get_alpha(data, estimator)

      # setup for individual estimates
      results <- data.frame(matrix(nrow = length(alpha), ncol = 7 + length(ICs)))
      names(results) <- c("alpha", "tuningparam", "tp", "DoF", "MSE",
                          "RMSFE_rec", "RMSFE_roll", ICs)
      results$alpha <- alpha

      # loop through alphas (~20 times)
      for (a in seq_along(alpha)){
        #print(paste("alpha =", alpha[a]))
      # 3) estimate
      estimates <- estimate(data, estimator, alpha[a])
      # estimate information
      results$tuningparam[a] <- estimates$tuning
      results$tp[a] <- estimates$tp
      results$DoF[a] <- estimates$DoF
      errsqrd <- (data$Y - estimates$My)^2
      results$MSE[a] <- sum(errsqrd)/data$T

      # # RMSFE TODO: In construction
      # results$RMSFE <- sqrt(mean((data$Y - estimates$My)^2))

      # 4) Information Criterias
      # GCV
      if ("GCV" %in% ICs){
        # Compute GCV
        numerator <- (1/data$T)*Norm(data$Y - estimates$My)^2
        denominator <- (1 - (1/data$T)*sum(diag(estimates$M)))^2
        results$GCV[a] <- numerator/denominator
      }
      # Mallows
      if ("Mallows" %in% ICs){
        # Compute Mallow's C_L
        residuals <- data$Y - data$X %*% estimates$delta
        sigma2hat <- sum((residuals - mean(residuals))^2) / data$T
        MIC <- Norm(data$Y - estimates$My)^2 / data$T + 2 * sigma2hat * sum(diag(estimates$M)) / data$T
        results$Mallows[a] <- MIC
        }
      # Leave One Out Cross Validation
      if ("LOOCV" %in% ICs){
        PSE <- c()
        for (i in 1:data$T){
          PSE[i] <- ((data$Y[i] - estimates$My[i])/1-estimates$M[i,i])^2
        }
        results$LOOCV[a] <- 1/data$T*sum(PSE)
      }
      # AIC
      if ("AIC" %in% ICs){
        Ahat <- as.matrix(estimates$eigen$vectors[,1:alpha[a]])
        Fhat <- t(Ahat) %*% data$X
        sqrdres <- (data$X - Ahat %*% Fhat)^2
        V <- sum(sqrdres) / (data$N * data$T)
        AIC <- log(V) + alpha[a] * 2 / data$T

        results$AIC[a] <- AIC
      }
      if ("AIC_sig" %in% ICs){
        Ahat <- as.matrix(estimates$eigen$vectors[,1:alpha[a]])
        Fhat <- t(Ahat) %*% data$X
        sqrdres <- (data$X - Ahat %*% Fhat)^2
        V <- sum(sqrdres) / (data$N * data$T)
        AIC_sig <- V + alpha[a] * 2 / data$T

        results$AIC_sig[a] <- AIC_sig
      }
      # BIC
      if ("BIC" %in% ICs){
        Ahat <- as.matrix(estimates$eigen$vectors[,1:alpha[a]])
        Fhat <- t(Ahat) %*% data$X
        sqrdres <- (data$X - Ahat %*% Fhat)^2
        V <- sum(sqrdres) / (data$N * data$T)
        BIC <- log(V) + alpha[a] * log(data$T) / data$T

        results$BIC[a] <- BIC
      }
      if ("BIC_sig" %in% ICs){
        Ahat <- as.matrix(estimates$eigen$vectors[,1:alpha[a]])
        Fhat <- t(Ahat) %*% data$X
        sqrdres <- (data$X - Ahat %*% Fhat)^2
        V <- sum(sqrdres) / (data$N * data$T)
        BIC_sig <- V + alpha[a] * log(data$T) / data$T

        results$BIC_sig[a] <- BIC_sig
      }
      # Bai-Ng
      if ("Bai_Ng_ic" %in% ICs){
        Ahat <- as.matrix(estimates$eigen$vectors[,1:alpha[a]])
        Fhat <- t(Ahat) %*% data$X
        sqrdres <- (data$X - Ahat %*% Fhat)^2
        V <- sum(sqrdres) / (data$N * data$T)
        ICp2 <- log(V) + estimates$tp * ((data$N + data$T) / (data$N *data$T)) * log(min(data$N, data$T))

        results$Bai_Ng_ic[a] <- ICp2
      }
      if ("Bai_Ng_pc" %in% ICs){
        Ahat <- as.matrix(estimates$eigen$vectors[,1:alpha[a]])
        Fhat <- t(Ahat) %*% data$X
        sqrdres <- (data$X - Ahat %*% Fhat)^2
        V <- sum(sqrdres) / data$N * data$T

        res <- data$X - Ahat %*% Fhat
        sighat <- sum((res - mean(res))^2) / (data$N * data$T)

        ICp2 <- V + estimates$tp * sighat * ((data$N + data$T) / (data$N *data$T)) * log(min(data$N, data$T))

        results$Bai_Ng_pc[a] <- ICp2
      }
    }


      # 5) Select optimal tuning parameters
      if (estimator %in% c("ridge", "LF")){
        results <- results %>% filter(Mallows == min(Mallows) | GCV == min(GCV))

        # Assign IC to optimal tuning parameter
        if (nrow(results) == 1){
          results <- rbind(results, results)
          results <- results %>% mutate(IC = c("Mallows", "GCV"))
        } else if (nrow(results) == 2){
          results <- results %>%
            mutate(IC = ifelse(Mallows == min(Mallows), "Mallows",
                         ifelse(GCV == min(GCV), "GCV", "fail")))

        }
      } else if (estimator == "PLS"){
        results <- results %>% filter(LOOCV == min(LOOCV)) %>%
          mutate(IC = "LOOCV")

      } else if (estimator == "PC"){
        results <- results %>% filter(Mallows == min(Mallows) |
                                      GCV == min(GCV) |
                                      AIC == min(AIC) |
                                      BIC == min(BIC) |
                                      AIC_sig == min(AIC_sig) |
                                      BIC_sig == min(BIC_sig) |
                                      Bai_Ng_ic == min(Bai_Ng_ic) |
                                      Bai_Ng_pc == min(Bai_Ng_pc)) %>%
          rename(Mallows_val = Mallows, GCV_val = GCV, AIC_val = AIC, AIC_sig_val = AIC_sig,
                 BIC_val = BIC, BIC_sig_val = BIC_sig, Bai_Ng_ic_val = Bai_Ng_ic, Bai_Ng_pc_val = Bai_Ng_pc)

        results <- results %>%
          mutate(Mallows = ifelse(Mallows_val == min(Mallows_val), 1, NA),
                 GCV = ifelse(GCV_val == min(GCV_val), 1, NA),
                 AIC = ifelse(AIC_val == min(AIC_val), 1, NA),
                 BIC = ifelse(BIC_val == min(BIC_val), 1, NA),
                 AIC_sig = ifelse(AIC_sig_val == min(AIC_sig_val), 1, NA),
                 BIC_sig = ifelse(BIC_sig_val == min(BIC_sig_val), 1, NA),
                 Bai_Ng_ic = ifelse(Bai_Ng_ic_val == min(Bai_Ng_ic_val), 1, NA),
                 Bai_Ng_pc = ifelse(Bai_Ng_pc_val == min(Bai_Ng_pc_val), 1, NA)) %>%
          pivot_longer(cols = "Mallows":"Bai_Ng_pc",
                       names_to = "IC", values_to = "drop") %>%
          drop_na(drop)

      }
      results$Estimator <- estimator
      results$DGP <- data$DGP
      results <- results %>% select(Estimator, DGP, IC, alpha, tuningparam,
                                    tp, DoF, MSE, RMSFE_rec, RMSFE_roll)

      # 6) Append optimal tuning parameters to output
      output <- rbind(output, results)
    }
  }
  return(output)
}


system.time(testPLS <- monte_carlo("large", "PLS", 1))
system.time(testRidge <- monte_carlo("large", "ridge", 1))
system.time(testLF <- monte_carlo("large", "LF", 1))
system.time(testPC <- monte_carlo("large", "PC", 1))






## =============================================================================
## ==  Run Simulations for Paper    ============================================
## =============================================================================

# ridge small
ridge_small <-  monte_carlo("small", "ridge", 10000)
save(ridge_small, file = paste0(simoutpath, "ridge_small.R"))

ridge_small_summary <- ridge_small %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
save(ridge_small_summary, file = paste0(simoutpath, "ridge_small_summary.R"))


# PLS large
PLS_large <-  monte_carlo("large", "PLS", 100)
save(PLS_large, file = paste0(simoutpath, "PLS_large.R"))

PLS_large_summary <- PLS_large %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
write_csv(PLS_large_summary, file = paste0(simoutpath, "PLS_large_summary.csv"))



# LF large
LF_large <-  monte_carlo("large", "LF", 100)
save(LF_large, file = paste0(simoutpath, "LF_large.R"))

LF_large_summary <- LF_large %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
write_csv(LF_large_summary, file = paste0(simoutpath, "LF_large_summary.csv"))



# 4) PC
# PC large slow
PC_large <-  monte_carlo("large", "PC", 100)
save(PC_large, file = paste0(simoutpath, "PC_large.R"))

PC_large_summary <- PC_large %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
save(PC_large_summary, file = paste0(simoutpath, "PC_large_summary.R"))



# PC small is fast
PC_small <-  monte_carlo("small", "PC", 10000)
save(PC_small, file = paste0(simoutpath, "PC_small.R"))

PC_small_summary <- PC_small %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
write_csv(PC_small_summary, file = paste0(simoutpath, "PC_small_summary.csv"))


# 3) LF
LF_large <-  monte_carlo("large", "LF", 100)
save(LF_large, file = paste0(simoutpath, "LF_large.R"))

LF_large_summary <- LF_large %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
write_csv(LF_large_summary, file = paste0(simoutpath, "LF_large_summary.csv"))

LF_small <-  monte_carlo("small", "LF", 10000)
save(LF_small, file = paste0(simoutpath, "LF_small.R"))

LF_small_summary <- LF_small %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
write_csv(LF_small_summary, file = paste0(simoutpath, "LF_small_summary.csv"))




# 2) Ridge
ridge_large <-  monte_carlo("large", "ridge", 100)
save(ridge_large, file = paste0(simoutpath, "ridge_large.R"))

ridge_large_summary <- ridge_large %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
save(ridge_large_summary, file = paste0(simoutpath, "ridge_large_summary.R"))

ridge_small <-  monte_carlo("small", "ridge", 1000)
save(ridge_small, file = paste0(simoutpath, "ridge_small.R"))

ridge_small_summary <- ridge_small %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
save(ridge_small_summary, file = paste0(simoutpath, "ridge_small_summary.R"))




# 1) PLS
PLS_large <-  monte_carlo("large", "PLS", 100)
save(PLS_large, file = paste0(simoutpath, "PLS_large.R"))

PLS_large_summary <- PLS_large %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
save(PLS_large_summary, file = paste0(simoutpath, "PLS_large_summary.R"))

PLS_small <-  monte_carlo("small", "PLS", 10000)
save(PLS_small, file = paste0(simoutpath, "PLS_small.R"))

PLS_small_summary <- PLS_small %>%
  group_by(Estimator, DGP, IC, tuningparam) %>%
  summarise(alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            tp_mean = mean(tp), tp_sd = sd(tp),
            DoF_mean = mean(DoF), DoF_sd = sd(DoF),
            MSE_mean = mean(MSE), MSE_sd = sd(MSE))
write_csv(PLS_small_summary, file = paste0(simoutpath, "PLS_small_summary.csv"))











