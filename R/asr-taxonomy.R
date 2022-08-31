attributeMatrix <- function(df, id_col, attr_col) {
    df[["Presence"]] <- 1
    select_cols <- c(id_col, unique(df[[attr_col]]))
    mat <- tidyr::pivot_wider(
        df, names_from = attr_col, values_from = "Presence", values_fill = 0
    )
    mat <- unique(mat[,select_cols])
    output <- as.data.frame(mat[,-1])
    rownames(output) <- mat[[id_col]]
    as.matrix(output)
}

estimateASR <- function(x, phy) {

    ## Ancestral charcter estimation
    ace_result <- ape::ace(
        x = x, phy = phy, type = "discrete",
        method = "ML", model = "ER", CI = TRUE,
        marginal = TRUE, kappa = 0
    )

    ## Extract likelihood as dataframe
    ace_states <- as.data.frame(round(ace_result$lik.anc, 3))

    ace_states
}


annotateInternalNodes <- function(ace, df, phy) {
    limit <- length(phy$tip.label) + 1
    all_nodes <- dplyr::arrange(df, node)
    internal_nodes <- all_nodes[limit:nrow(all_nodes),]
    new_attribute_data <- cbind(internal_nodes, ace)
    output <- new_attribute_data[!is.na(new_attribute_data$taxid),]
    output$Evidence <- "ASR"
    output$Confidence_interval <- dplyr::case_when(
        output$`1` == 1 ~ "always",
        output$`1` < 1 & output$`1` >= .9 ~ "usually",
        output$`1` < .9 & output$`1` >= .8 ~ "normally",
        output$`1` < .8 & output$`1`>= .7 ~ "often",
        output$`1` < .7 & output$`1` >= .5 ~ "sometimes",
        output$`1` < .5 & output$`1` >= .3 ~ "occasionally",
        output$`1` < .3 & output$`1` >= .1 ~ "seldom",
        output$`1` < .1 & output$`1` >= .05 ~ "rarely",
        output$`1` < .05 & output$`1` >= 0 ~ "never"

    )
    output$Taxon_name <- gsub("_", " ", output$sci_name)
    output$NCBI_ID <- output$taxid
    output$Rank <- output$rank
    output$Attribute_source <- TRUE
    select_cols <- c(
        "NCBI_ID", "Taxon_name", "Evidence",
        "Confidence_interval", "Rank"
    )
    output[, select_cols]
}


createEntriesASR <- function(x, phy, df) {

    valid_ranks <- c(
        "superkingdom", "phylum", "class", "order", "family", "genus",
        "species", "strain"
    )

    output <- vector("list", length(x))
    names(output) <- names(x)
    for (i in seq_along(output)) {
        output[[i]] <- x[[i]] %>%
            estimateASR(phy = phy) %>%
            annotateInternalNodes(df = df, phy = phy)
    }
    output %>%
        dplyr::bind_rows(.id = "Attribute") %>%
        dplyr::filter(
            Confidence_interval != "never",
            Rank %in% valid_ranks) %>%
        dplyr::relocate(
            NCBI_ID, Taxon_name, Attribute, Evidence,
            Confidence_interval, Rank
        )
}

appendASR <- function(df1, df2) {
    df1_taxids <- df1$NCBI_ID
    df2 <- df2[!df2$NCBI_ID %in% df1_taxids, ]
    dplyr::bind_rows(df1, df2)
}

makeSignatures <- function(
        df,
        tax.id.type = "ncbi",
        tax.level = "mixed",
        min.size = 1,
        Evidence = c("EXP", "ASR"),
        Confidence_interval = 1
) {
    valid_ranks <- c(
        "superkingdom", "phylum", "class", "order", "family", "genus",
        "species", "strain"
    )

    if (tax.level == "mixed") {
        tax.level <- valid_ranks
    }

    if (tax.id.type == "ncbi") {
        tax.id.type <- "NCBI_ID"
    } else if (tax.id.type == "taxname") {
        tax.id.type <-  "Taxon_name"
    }

    ci_vals <- c(
        "always", "usually", "normally", "often", "sometimes", "occasionally",
        "seldom", "rarely"
    )

    if (Confidence_interval == 1) {
        Confidence_interval <- ci_vals[1]
    } else if (Confidence_interval < 1 & Confidence_interval >= .9) {
        Confidence_interval <- ci_vals[1:2]
    } else if (Confidence_interval < .9 & Confidence_interval >= .8) {
        Confidence_interval <- ci_vals[1:3]
    } else if (Confidence_interval < .8 & Confidence_interval >= .7) {
        Confidence_interval <- ci_vals[1:4]
    } else if (Confidence_interval < .7 & Confidence_interval >= .5) {
        Confidence_interval <- ci_vals[1:5]
    } else if (Confidence_interval < .5 & Confidence_interval >= .3) {
        Confidence_interval <- ci_vals[1:6]
    } else if (Confidence_interval < .3 & Confidence_interval >= .1) {
        Confidence_interval <- ci_vals[1:7]
    } else if (Confidence_interval < .1 & Confidence_interval >= .05) {
        Confidence_interval <- ci_vals[1:8]
    } else if (Confidence_interval < .05 & Confidence_interval >= 0) {
        Confidence_interval <- ci_vals[1:8]
    }

    df <- df %>%
        dplyr::filter(
            .data[["Rank"]] %in% .env[["tax.level"]],
            .data[["Evidence"]] %in% .env[["Evidence"]],
            .data[["Confidence_interval"]] %in% .env[["Confidence_interval"]]
        )

    select_cols <- c(tax.id.type, "Attribute")
    output <- df[,select_cols]
    output <- split(output, factor(output[["Attribute"]]))
    output <- lapply(output, function(x) {
        values <- x[[1]]
        names <- x[[2]]
        names(values) <- names
        unique(values)
    })
    output[vapply(output, function(x)  length(x) >= min.size, logical(1))]
}
