
x <- m3
x[x > 0] <- 1
x <- x[,sort(colnames(x))]
rownames(x) <- NULL
y <- as.matrix(unique(as.data.frame(x)))
rownames(y) <- NULL
output <- vector('character', nrow(y))
for (i in seq_along(output)) {
    pos <- which(y[i,] == 1)
    attr_names <- colnames(y)[pos]
    output[[i]] <- paste0(attr_names, collapse = ':')
}
rownames(y) <- output

trait_mat <- m3[, sort(colnames(m3))]
trait_mat[trait_mat > 0] <- 1
trait_v <- vector('character', nrow(trait_mat))
for (i in seq_along(trait_v)) {
    pos <- which(trait_mat[i,] == 1)
    attr_names <- colnames(trait_mat)[pos]
    trait_v[i] <- paste0(attr_names, collapse = ':')
}

myTraits <- matrix(
    data = trait_v, ncol = 1, dimnames = list(
        rownames = rownames(trait_mat),
        colnames = 'trait'
    )
)

myData <- phyDat(
    myTraits, type = "USER",
    levels = rownames(y),
    contrast = y
)

myData

system.time(
    anc_pars <- ancestral.pars(mpa_tree, myData, type="MPR")
)


anc_pars_data <- vector('list', length(anc_pars))
for (i in seq_along(anc_pars_data)) {
    anc_pars_data[[i]] <- as.data.frame(unname(anc_pars[[i]]))
}


row_names <- c(mpa_tree$tip.label, mpa_tree$node.label)
col_names <-
final <- bind_rows(anc_pars_data)
rownames(final) <- row_names
colnames(final) <- colnames(y)

hist(as.matrix(final))
