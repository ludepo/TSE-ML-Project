# TSE-ML-Project

This repository contains the code written for the final project of the class "Machine Learning for Economics" during my masters program at the Toulouse School of Economics. The project was done in cooperation with Nikita Marini and includes a replication of the simulations presented in Carrasco & Rossi (2016) as well as a practical application using the Stock & Watson (2002) dataset on macroeconomic variables. 

## Overview

The Carrasco & Rossi (2016) paper focused on in-sample prediction and out-of-sample forecasting in high dimensional models with many exogenous regressors. Four dimension reduction devices are used: Principal Components, Ridge Regression, Landweber Fridman (LF), and Partial Least Squares (PLS). Each involves a regularization or tuning parameter that is selected through generalized cross validation (GCV), Mallows Cp, Leave one out cross validation (LOOCV), AIC, BIC or the Bai & Ng (2002) IC and PC criterion. Following Carrasco and Rossi we evaluate these estimators in a monte carlo simulation framework with 6 different data generating processes (DGPs) and a subsequent application on the macroeconomic dataset constructed by Stock & Watson (2002).

## Structure
The 'CODE' folder contains the following scripts: 
* _Master:_ master script to replicate the Carrasco & Rossi (2016) simulation, calls the 'Estimators' and 'Simulation' script
  + _Simulation:_ contains function to create simulated data according to different DGPs
  + _Estimators:_ contains fuction for determining the alpha range as well as function with estimators
* _SW_data:_ script to fetch the Stock & Watson (2002) dataset used for the empirical application
* _Application:_ scirpt for empirical application
