# Objective : Apply estimators to Stock and Watson Data
# Created by: Luca Poll
# Created on: 4/17/2021


library(tidyverse)

# define paths
#setwd("..")
scripts <- paste0(getwd(), "/CODE/")
datapath <- paste0(getwd(), "/DATA/")
simoutpath <- paste0(getwd(), "/OUTPUT/SIMULATION/")

# load functions
source(paste0(scripts, "Simulation.R"))
source(paste0(scripts, "Estimators.R"))



# application functions
appl_alpha <- function(data, estimator){

    # get alpha range for ridge estimator
    if (estimator == "ridge"){
      alpha <- seq(0.000000001, 3.5, length = 100)
    }

    # get alpha range for Landweber Fridman estimator
    if (estimator == "LF"){
      alpha <- seq(0.00000000001, 0.01, length = 100)
    }

	if (estimator == "PC"){
		alpha <- 1:data$rmax
	}

    if (estimator == "PLS"){

		alpha <- seq(1, 35)
    }

    return(alpha)

}


find_min <- function(estimator, data){

  # setup
  output <- data.frame(Estimator = character(), IC = character(),
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
    ICs <- c("GCV", "Mallows", "AIC", "BIC", "AIC_sig", "BIC_sig", "Bai_Ng_ic")
  }


  # 2) get alpha range
  alpha <- appl_alpha(data, estimator)

  # setup for individual estimates
  results <- data.frame(matrix(nrow = length(alpha), ncol = 7 + length(ICs)))
  names(results) <- c("alpha", "tuningparam", "tp", "DoF", "MSE",
                          "RMSFE_rec", "RMSFE_roll", ICs)
  results$alpha <- alpha

  # loop through alphas (~20 times)
  for (a in seq_along(alpha)){
      print(paste("alpha =", alpha[a]))
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

  results$Estimator <- estimator
  results <- results %>% select(Estimator, all_of(ICs), alpha, tuningparam,
                                    tp, DoF, MSE, RMSFE_rec, RMSFE_roll)

  # 6) Append optimal tuning parameters to output
  output <- rbind(output, results)

  # 7) Plot IC for values of alpha
  plot_data <- output %>% select(all_of(ICs), alpha, DoF) %>%
      pivot_longer(!c(alpha, DoF), names_to = "IC", values_to = "criterion")

  if (estimator == "PLS"){
    plot <- ggplot(plot_data, aes(x=alpha, y=criterion, group=IC)) +
      geom_line(aes(color=IC), size=1) + theme_minimal() +
      geom_point(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion)),
                 aes(x=alpha, y=criterion, color=IC), size=3) +
      geom_text(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion)),
                aes(x=alpha, y=criterion, color=IC,
                    label = paste("criterion =", round(criterion, 6),
                                  "\n r =", alpha)),
                vjust=-0.5, hjust=-.1,size=4) +
      scale_color_brewer(palette = "Dark2")
    plot2 <- NULL
  }

  if (estimator %in% c("ridge", "LF")){
    plot <- ggplot(plot_data, aes(x=alpha, y=criterion, group=IC)) +
    geom_line(aes(color=IC), size=1) + theme_minimal() +
    geom_point(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion)),
               aes(x=alpha, y=criterion, color=IC), size=3) +
    geom_text(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion)),
              aes(x=alpha, y=criterion, color=IC,
                  label = paste(IC, "criterion",
                                "\n DoF =", round(DoF))),
              vjust=-0.5, hjust=-.1,size=4) +
    scale_color_brewer(palette = "Dark2")
    plot2 <- NULL
  }
  if (estimator == "PC"){
    plot <- ggplot(plot_data, aes(x=alpha, y=criterion, group=IC)) +
      geom_line(aes(color=IC), size=1) + theme_minimal() +
      geom_point(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion)),
                 aes(x=alpha, y=criterion, color=IC), size=3) +
      geom_text(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion)),
                aes(x=alpha, y=criterion, color=IC,
                    label = paste(IC, "criterion",
                                  "\n r =", alpha)),
                vjust=-0.5, hjust=-.1,size=4) +
      xlim(0, max(alpha)+max(alpha)/6) + scale_color_brewer(palette = "Dark2")

    plot2 <- plot_data %>% filter(IC %in% c("GCV", "Mallows")) %>%
      ggplot(aes(x=alpha, y=criterion, group=IC)) +
      geom_line(aes(color=IC), size=1) + theme_minimal() +
      geom_point(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion),
                                                              IC %in% c("GCV", "Mallows")),
                 aes(x=alpha, y=criterion, color=IC), size=3) +
      geom_text(data = plot_data %>% group_by(IC) %>% filter(criterion == min(criterion),
                                                             IC %in% c("GCV", "Mallows")),
                aes(x=alpha, y=criterion, color=IC,
                    label = paste(IC, "criterion",
                                  "\n r =", alpha)),
                vjust=-0.5, hjust=-.1,size=4) +
      xlim(0, max(alpha)+max(alpha)/6) + scale_color_brewer(palette = "Dark2")
  }

  return(list(output = output, plot = plot, plot2 = plot2))
}


## -----------------------------------------------------------------------------
## --- Load Data and run model     ---------------------------------------------
## -----------------------------------------------------------------------------

load(paste0(datapath, "SW_data.R"))

ridge <- find_min("ridge", data)
ridge$plot

LF <- find_min("LF", data)
LF$plot

PLS <- find_min("PLS", data)
PLS$plot

PC <- find_min("PC", data)
PC$plot
PC$plot2



# paste results together in table
rtable <- ridge$output %>%
  filter(Mallows == min(Mallows) | GCV == min(GCV)) %>%
  mutate(IC = ifelse(Mallows == min(Mallows), "Mallows",
                     ifelse(GCV == min(GCV), "GCV", "fail")),
         DoF = round(DoF)) %>%
  select(Estimator, IC, tuningparam, tp, DoF, MSE)

lftable <- LF$output %>%
  filter(Mallows == min(Mallows) | GCV == min(GCV)) %>%
  mutate(IC = ifelse(Mallows == min(Mallows), "Mallows",
                     ifelse(GCV == min(GCV), "GCV", "fail")),
         DoF = round(DoF)) %>%
  select(Estimator, IC, tuningparam, tp, DoF, MSE)

plstable <- PLS$output %>% filter(LOOCV == min(LOOCV)) %>%
          mutate(IC = "LOOCV") %>%
  select(Estimator, IC, tuningparam, tp, DoF, MSE)

pctable1 <- PC$output %>% filter(Mallows == min(Mallows) | GCV == min(GCV) |
                                  AIC == min(AIC) | BIC == min(BIC) |
                                  AIC_sig == min(AIC_sig) | BIC_sig == min(BIC_sig) |
                                  Bai_Ng_ic == min(Bai_Ng_ic)) %>%
          rename(Mallows_val = Mallows, GCV_val = GCV, AIC_val = AIC, AIC_sig_val = AIC_sig,
                 BIC_val = BIC, BIC_sig_val = BIC_sig, Bai_Ng_ic_val = Bai_Ng_ic)
pctable <- pctable1 %>%
  mutate(Mallows = ifelse(Mallows_val == min(Mallows_val), 1, NA),
         GCV = ifelse(GCV_val == min(GCV_val), 1, NA),
         AIC = ifelse(AIC_val == min(AIC_val), 1, NA),
         BIC = ifelse(BIC_val == min(BIC_val), 1, NA),
         AIC_sig = ifelse(AIC_sig_val == min(AIC_sig_val), 1, NA),
         BIC_sig = ifelse(BIC_sig_val == min(BIC_sig_val), 1, NA),
         Bai_Ng_ic = ifelse(Bai_Ng_ic_val == min(Bai_Ng_ic_val), 1, NA)) %>%
  pivot_longer(cols = "Mallows":"Bai_Ng_ic", names_to = "IC", values_to = "drop") %>%
  drop_na(drop) %>% mutate(DoF = NA) %>%
  select(Estimator, IC, tuningparam, tp, DoF, MSE)

results <- rbind(rtable, lftable) %>%
  rbind(plstable) %>%
  rbind(pctable)

library(xtable)

xtable(results)

library(factoextra)
res.pca <- prcomp(data$X, scale = TRUE)
fviz_eig(res.pca)


