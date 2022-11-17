library(bugphyzz)
aer <- physiologies('aerophilicity')[[1]]
data_backup <- propagateUpstream(aer)

data <- data_backup
valid_ranks <- .validRanks()
valid_ranks <- valid_ranks[1:which(valid_ranks == 'family')]

for (i in seq_along(valid_ranks)) {

    current_rank <- valid_ranks[i]

    if (current_rank %in% names(data)) {
        message('Current rank: ', current_rank)
        next_pos <- i + 1
        if (next_pos > length(valid_ranks))
            break()
        next_rank <- valid_ranks[next_pos]
        message('Getting next rank: ', next_rank)
        parent_scores <- getParentScores(data[[current_rank]])
        parent_scores <- parent_scores |>
            dplyr::filter(.data$Rank == next_rank)

        if (next_rank %in% names(data)) {

            data[[next_rank]] <-
                .replaceParents(data[[next_rank]], parent_scores)

        } else {
            data[[next_rank]] <- parent_scores
        }

    } else {
        next()
    }
}


res <- propagateUpstream(aer, 'family')

lapply(data, function(x) unique(x$Evidence))

