# Title

Quasi-maximum Likelihood Inference for Linear Double Autoregressive Models



# Description of this file

- Dataset
  - **BTCweekly**: BTCweekly: The named 'BTCweekly' dataset contains weekly closing prices of the Bitcoin (BTC) from July 18, 2010 to August 16, 2020. This dataset has two variables: Date and closing price, which represent the date and BTC respectively. 
- Programs
  - **rcpp_LDAR.cpp**: Basic functions in form of Rcpp code.
  - **LDAR.R**: R functions to perform the proposed estimation and inference as well as to support 'rcpp_LDAR.cpp'.
  - **KSD.R**: R functions to conduct the Kernelized Stein Discrepancy (KSD) test.
  - **DWQRE.R**: R functions to calculate the doubly weighted quantile regression estimator (DWQRE).
  
- **Rcode_for_LDAR.html**: This file contains all outputs for real data analysis of Liu et al. (2022+). Specifically, (1) Model estimations of the LDAR model based on EQMLE, GQMLE and DWQRE are presented in each subsection; (2) Based on the residuals of E-QMLE, we use Q-Q plot and KSD test to check the error distribution, and conduct portmanteau test to check the adequacy of fitted models; (3) The AREs are provided to compare the efficiency of three estimation methods. (4) Rolling forecasting and backtests are presented.

- **Rcode_for_LDAR.Rmd**: Source code for **Rcode_for_LDAR.html**.



# Reference

Liu, H., Tan, S., and Zhu, Q.(2022+) Quasi-maximum Likelihood Inference for Linear Double Autoregressive Models.



# Contact information
Name: Tan Songhua

Email: tansonghua@163.sufe.edu.cn 
