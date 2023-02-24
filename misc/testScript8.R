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
