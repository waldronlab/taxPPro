
library(taxPPro)
library(data.tree)
library(ape)

data("tree_list")
data_tree <- as.Node(tree_list)
Prune(data_tree, function(x) !grepl("t__", x$name))
phylo <- as.phylo.Node(x = data_tree, heightAttribute = NULL)
phyloResolvedMulti <- multi2di(phylo)
phyloCollapsedSingles <- collapse.singles(phyloResolvedMulti, root.edge = TRUE)
phyloFname <- file.path('inst', 'extdata', 'ncbi_sp_tree_resolved.newick')
write.tree(phy = phyloCollapsedSingles, file = phyloFname)
