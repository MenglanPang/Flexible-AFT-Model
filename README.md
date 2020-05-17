# Flexible AFT Model
R code for implementing the method in article “Pang M, Platt R, Schuster T, Abrahamowicz M. Flexible Modeling of Time-dependent and Non-linear covariate effects in Accelerated FailureTime model. Under review at *Statistical Method in Medical Research* 2020".

## Content
## Description
A flexible extension of the accelerated failure time model that models:
- baseline hazard
- TD effects of continuous and categorical variables
- NL effects of continuous variables

The code has been written using R with the following version information:<br/>
- R version 3.3.1 (2016-06-21)<br/> 
- Platform x86_64-w64-mingw32/x64 (64-bit)<br/> 
- Using R packages:<br/> 
  - survival version 3.1-12
  - splines version v3.6.2
  
Code to implement the flexible AFT model:
#### `FlexAFT.R`
This program includes all necessary functions to provide estimates of
- Hazard function and survival curve conditional on an arbitrary covariate pattern
- TD effects
- NL effects
The program is called by the program `Example.R`. 

Code to run the flexible AFT model:
#### `Example.R`
The program read the data in `sepsis.csv`, and generates the results save in  `sepsis_fit.rda`
 
For questions or comments about the code please contact Menglan Pang (menglan.pang at mail.mcgill.ca).
