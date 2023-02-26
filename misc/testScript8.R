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
ph$Score <- sample(c(0.5, 1), length(ph$Score), replace = TRUE)

## range
w <- phys$width |> # numeric range
    prepareData() |>
    prepareData2()

data('tree_list')
tree <- data.tree::as.Node(tree_list)

x <- addAttributes(tree, aer)
x$Do(asrUpstreamLogical, traversal = 'post-order')
x$Do(inhDownstreamLogical, traversal = 'pre-order')

y <- addAttributes(tree, ph)
y$Do(asrUpstreamNumeric, traversal = 'post-order')
y$Do(inhDownstreamNumeric, traversal = 'pre-order')

z <- addAttributes(tree, w)

printDataTreeAttributes(z, limit = 1000) |> View()
