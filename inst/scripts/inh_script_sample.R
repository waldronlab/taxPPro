
kingdom <- c(
    aerobic = 0.7, anaerobic = 0, fac = 0.3, aerotolerant = 0
)

ranks <- c(
    'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'
)

output <- vector('list', length(ranks))
names(output) <- ranks
for (i in seq_along(output)) {
    if (i == 1) {
        output[[i]] <- taxPPro:::myFun(kingdom, adjF = 0.1)
    } else {
        output[[i]] <- taxPPro:::myFun(output[[i - 1]], adjF = 0.1)
    }
}
output <- c(list(kingdom = kingdom), output)
df <- purrr::map(output , ~ as.data.frame(matrix(round(.x, 2), nrow  = 1))) |>
    dplyr::bind_rows(.id = 'rank')
colnames(df) <- c('rank', 'aer', 'ana', 'fac', 'aet')
