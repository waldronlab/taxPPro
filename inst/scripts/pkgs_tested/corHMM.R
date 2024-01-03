library(corHMM)

data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
# # one way to get the parameters from your corHMM object in the correct order
p <- sapply(1:max(MK_3state$index.mat, na.rm = TRUE), function(x)
    na.omit(c(MK_3state$solution))[na.omit(c(MK_3state$index.mat) == x)][1])
# using custom params
states_1 <- ancRECON(phy = phy, data = MK_3state$data, p = p, method = "marginal",
                     rate.cat <- MK_3state$rate.cat, ntraits = NULL, rate.mat = MK_3state$index.mat,
                     root.p = MK_3state$root.p)
