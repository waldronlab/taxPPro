library(bugphyzz)

phys <- physiologies(
    keyword = c('aerophilicity', 'optimal ph', 'width'),
    remove_false = TRUE, full_source = FALSE
)

## logical
aer <- phys$aerophilicity |>
    prepareData() |>
    prepareData2()

## numeric
ph <- phys$`optimal ph` |> # numeric
    prepareData() |>
    prepareData2()

## range
w <- phys$width |> # numeric range
    prepareData() |>
    prepareData2()

data('tree_list')
tree <- data.tree::as.Node(tree_list)
x <- addAttributes(tree, aer)
y <- addAttributes(tree, ph)
z <- addAttributes(tree, w)

View(printDataTreeAttributes(y, limit = 100))
