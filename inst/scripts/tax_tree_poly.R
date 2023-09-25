phys_data_ready



mat <- phys_data_ready |>
    select(NCBI_ID, Attribute, Score) |>
    filter(!is.na(Attribute)) |>
    # complete(NCB_ID, Attribute, fill = list(Score = 0)) |>
    pivot_wider(
        names_from = Attribute, values_from = 'Score', values_fill = 0
    ) |>
    tibble::column_to_rownames(var = 'NCBI_ID') |>
    as.matrix()




data('tree_list')
data_tree <- as.Node(tree_list)
phylo_tree <- as.phylo.Node(data_tree, heightAttribute = NULL)

phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))


myIndex <- which(!phylo_tree$tip.label %in% rownames(mat))
no_annotated_tips <- phylo_tree$tip.label[myIndex]



mat2 <- matrix(
    data = rep(rep(1/ncol(mat), ncol(mat)), length(no_annotated_tips)),
    nrow = length((no_annotated_tips)),
    byrow = TRUE,
    dimnames = list(rownames = no_annotated_tips, colnames = colnames(mat))
)

mat3 <- rbind(mat, mat2)
mat3 <- mat3[phylo_tree$tip.label,]
dim(mat3)

myFit <- fitMk(tree = phylo_tree, x = mat3, model = 'ER')
myAsr <- ancr(object = myFit)



