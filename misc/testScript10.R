
library(bugphyzz)
library(purrr)
library(dplyr)

phys_names <- c('aerophilicity', 'growth temperature')
phys <- physiologies(phys_names, remove_false = TRUE, full_source = TRUE)
aer <- phys$aerophilicity
gt <- phys$`growth temperature`
a <- getDataReadyForPropagation(aer)
g <- getDataReadyForPropagation(gt)
a_ag <- getAgreements(a)
g_ag <- getAgreements(g)
x <- resolveAgreements(a)
y <- resolveAgreements(g)
j <- resolveConflicts(x)
k <- resolveConflicts(y)
