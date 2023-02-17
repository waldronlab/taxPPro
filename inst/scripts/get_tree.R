library(taxPPro)
library(data.tree)
tbl <- getNCBI()
tree <- as.Node(tbl[, 'pathString', drop = FALSE], pathDelimiter = '|||')
tree$Do(addStrains)

print(tree, limit = 1000) |> View()

