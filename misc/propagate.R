
library(bugphyzz)
library(taxPPro)
# library(dplyr)
# library(purrr)
# library(tidyr)

# aer <- tibble::as_tibble(physiologies('aerophilicity')[[1]])
# aer_extended <- propagate(aer)

# ph <- physiologies('optimal ph')[[1]]
# ph_extended <- propagate(ph)

# filtered_ph <- ph |>
#     filter_dataset_for_propagation() |>
#     remove_taxa_duplicates() |>
#     ci_to_scores() |>
#     distinct()
#
# x3 <- calcParentScore(filtered_ph)
#
# x2 <- ph |>
#     upstream()

# k -----------------------------------------------------------------------

calcScore <- function(df, wt = TRUE) {

    attr <- val <- NULL
    colnames(df) <- c('attr', 'val')

    if (wt) {
        df |> mutate(
            total = sum(val),
            prop = val / total
        ) |>
            dplyr::count(attr, wt = prop, name = 'Score') |>
            dplyr::mutate(Score = round(Score, 1)) |>
            dplyr::filter(Score >= 0.5)
    }
}


x <- c(A = 1, A = 1, A = 1, B = 0.5)

df <- data.frame(
    x = c(rep('sp1', 4), rep('sp2', 4)),
    y = c(rep('A',3), 'B', rep('A', 3), 'B'),
    z = c(rep(1, 3), 0.5, rep(1, 3), 0.5)
)

df2 <- df |>
    group_by(x) |>
    tidyr::nest()

df2 |>
    mutate(data2 = purrr::map(data, calcScore)) |>
    select(-data) |>
    tidyr::unnest(data2)



