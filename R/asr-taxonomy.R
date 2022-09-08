
#' Estimate ASR
#'
#' \code{extimateASR} performs ASR.
#'
#' @param x A a named character vector. The values are the observed character
#' states and the names are the tip labels. It must be of the same length as
#' the number of tip labels in the phylogenetic tree.
#' @param phy A phylo object. With tip labels and branch lenghts.
#'
#' @return A dataframe.
#' @export
#'
estimateASR <- function(x, phy) {

    ## Ancestral character estimation
    ace_result <- ape::ace(
        x = x, phy = phy, type = "discrete",
        method = "ML", model = "ER", CI = TRUE,
        marginal = TRUE, kappa = 0
    )

    ## Extract likelihood as dataframe
    ace_states <- as.data.frame(round(ace_result$lik.anc, 3))

    ace_states
}


#' Annotate internal nodes
#'
#' \code{annotateInternalNodes} annotates the internal nodes of the output of
#' the `estimateASR` function.
#'
#' @param ace Output of the estimateASR function.
#' @param df A data frame.
#' @param phy A phylo object.
#'
#' @return A data frame
#' @export
#'
annotateInternalNodes <- function(ace, df, phy) {
    limit <- length(phy$tip.label) + 1
    all_nodes <- dplyr::arrange(df, node)
    internal_nodes <- all_nodes[limit:nrow(all_nodes),]
    new_attribute_data <- cbind(internal_nodes, ace)
    output <- new_attribute_data[!is.na(new_attribute_data$taxid),]
    output$Evidence <- "asr"
    output$Frequency <- dplyr::case_when(
        output$`1` == 1 ~ "always",
        output$`1` >= .7 & output$`1` < 1 ~ "usually",
        output$`1` >= .4 & output$`1` < .7 ~ "sometimes",
        output$`1` >= .1 & output$`1` < .4 ~ "rarely",
        output$`1` < .1 ~ "never"
    )
    output$Taxon_name <- gsub("_", " ", output$sci_name)
    output$NCBI_ID <- output$taxid
    output$Rank <- output$rank
    output$Attribute_source <- TRUE
    select_cols <- c(
        "NCBI_ID", "Taxon_name", "Evidence",
        "Frequency", "Rank"
    )
    output[, select_cols]
}

#' Create entries for ASR
#'
#' @param x Character states.
#' @param phy A phylo object.
#' @param df Tree metadata (a data frame).
#'
#' @return A data frame
#' @export
#'
createEntriesASR <- function(x, phy, df) {

    valid_ranks <- c(
        "superkingdom", "phylum", "class", "order", "family", "genus",
        "species", "strain"
    )

    output <- vector("list", length(x))
    names(output) <- names(x)
    for (i in seq_along(output)) {
        output[[i]] <- x[[i]] |>
            {\(y) estimateASR(x = y, phy = phy)}() |>
            {\(y) annotateInternalNodes(ace = y, df = df, phy = phy)}()
    }

    output %>%
        dplyr::bind_rows(.id = "Attribute") %>%
        dplyr::filter(
            Frequency != "never",
            Rank %in% valid_ranks
        ) %>%
        dplyr::relocate(
            NCBI_ID, Taxon_name, Attribute, Evidence,
            Frequency, Rank
        )
}

# appendASR <- function(df1, df2) {
#     df1_taxids <- df1$NCBI_ID
#     df2 <- df2[!df2$NCBI_ID %in% df1_taxids, ]
#     dplyr::bind_rows(df1, df2)
# }
