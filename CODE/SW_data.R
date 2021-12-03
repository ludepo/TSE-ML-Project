# Objective : Get Stock and Watson data
# Created by: Luca Poll
# Created on: 4/18/2021

library(tidyverse)
# install.packages("remotes")
# remotes::install_github("johannestang/macrods")
library(macrods)
library(quantmod)
library(VIM)
library(mice)
library(zoo)
library(psych)
library(reshape2)
library(Hmisc)
library(mice)


### Covariates

# Get Stock and Watson transformed data
SW <- getmacrodata(SW2002, enddate = 2000)

# dataframe
X <- cbind(data.frame(as.matrix(SW)), data.frame(date = as.Date(as.yearmon(time(SW))))) %>% relocate(date)

# missing
sort(sapply(X, function(x) sum(is.na(x))))

clean <- X %>% select(-c(FSTE, FSTM, FTMD, FSTB, FTB, CCINRV, CCINT, CCINV, FCLNBF, FWAFIT))

mice_plot <- aggr(clean, col=c("blue", "red"), combined = T, numbers=TRUE, labels = T, cex.numbers=0.8, cex.axis=.6)
matrixplot(clean)

# drop last year
cleaned <- clean %>% filter(date < as.Date("1999-01-01"))

# impute using predictive mean matching
impu <- mice(cleaned, m=1, maxit = 50, method='pmm')
compl <- complete(impu)

matrixplot(compl)
psych::describe(compl)

# standardize
compl_nor <- data.frame(scale(compl[,2:ncol(compl)]))
compl_nor <- cbind(compl_nor, compl[,1]) %>% rename(date = "compl[, 1]") %>% relocate(date)
psych::describe(compl_nor)

# correlation
cormat <- round(cor(compl_nor[,2:ncol(compl_nor)]),2)
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + theme_classic()



### Dependent variable

# Get S&P 500 adjusted prices and aggregate to quarterly
SPX <- data.frame(getSymbols("^GSPC",auto.assign = FALSE, from = "1958-12-01"))
SPX$date <- rownames(SPX)
sp <- SPX %>% mutate(month = str_sub(date, 6,7), year =  str_sub(date, 1,4)) %>% group_by(year, month) %>%
  summarise(price = mean(GSPC.Adjusted)) %>% mutate(date = paste(year, month, "01", sep = "-")) %>% ungroup() %>%
  mutate(logret  = log(price) - log(lag(price))) %>% select(date, logret) %>% drop_na() %>% mutate(date = as.Date(date))



### Combine
compl_nor <- compl_nor %>% left_join(sp, by = "date") %>% relocate(date, logret)

covariates <- as.matrix(compl_nor[,3:ncol(compl_nor)])
dep <- as.vector(compl_nor$logret)

data <- list(X = covariates, Y = dep, rmax = 100, N = dim(covariates)[2], T = dim(covariates)[1])

save(data, file = "C:/Users/Luca Poll/Google Drive/TSE/3rd Semester/Machine Learning/Project/TSE-ML-Project/DATA/SW_data.R")






