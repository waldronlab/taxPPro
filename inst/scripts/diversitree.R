library(diversitree)

## It seems that is based on the same models
## as ape::ace (Pagel 1994)

## Check out this link:
## https://www.zoology.ubc.ca/~fitzjohn/diversitree.docs/asr-bisse.html

pars <- c(0.1, 0.2, 0.03, 0.06, 0.01, 0.02)
set.seed(3)
phy <- trees(pars, "bisse", max.taxa = 50, max.t = Inf, x0 = 0)[[1]]

lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars, method = "subplex")
st <- asr.marginal(lik, coef(fit))
nodelabels(thermo = t(st), piecol = 1:2, cex = 0.5)
