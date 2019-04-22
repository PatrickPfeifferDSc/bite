# bite
Public repository for the R-package 'bite', abbreviation for Bayesian Inference on Treatment Effects

This R-package provides methods and functionalities of treatment effect estimation on outcomes measured in 
panel structure. The modelling strategies were originally proposed in 
Jacobi, Wagner, Frühwirth-Schnatter (2016),
"Bayesian treatment effects models with variable selection for panel data outcomes 
with an application to earnings effects of maternity leave", Journal of Econometrics, 193, 234-250.
There are little differences in methodology to the paper, except for what were necessary changes to implement
them in R.

For using this package you may refer to devtools:
> install.packages("devtools") \
> devtools::install_github("PatrickPfeifferDsc/bite")
Or CRAN:
install.packages("bite")
which should make little difference, except changes will be available on Github first.

If you feel something is wrong, or could be improved, comments are always appreciated.
Remember this is work in progress, still I will be happy to consider any comment.
