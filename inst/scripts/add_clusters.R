library(data.tree)
library(taxPPro)
library(purrr)
library(dplyr)

data('tree_list')
ncbiTree <- as.Node(tree_list)
g <- ncbiTree$Get(
    function(node) node,
    filterFun = function(node) grepl('g__', node$name)
)

cl <- map(g, ~ {
    .x |> 
        as.list() |> 
        unlist() |> 
        names() |> 
        strsplit('\\.') |> 
        unlist() |> 
        unique() |> 
        {\(y) sort(y[grep('[gst]__', y)])}() |> 
        paste0(collapse = '+')
}) |> 
    flatten_chr()

cl <- map2_chr(.x = cl, .y = names(g), ~ {paste0(.y, '+', .x)})

tips <- getDistantTips()

tip_data <- ltp()$tip_data

o <- map_chr(tip_data$NCBI_ID, ~ {
    if (is.na(.x)) 
        return(NA)
    regex <- paste0('\\b', .x, '\\b')
    res <- grep(regex, cl, value = TRUE)
    if (!length(res))
        return(NA)
    return(res)
})

tip_data$cluster <- o

myData <- tip_data[,c('tip_label', 'cluster')]

d <- getDistantTips()
d <- left_join(d, myData, by = c('tip1' = 'tip_label'))

write.table(
    x = d, file = '~/Projects/taxPPro/inst/extdata/longest_distance_between_tips.tsv',
    sep = '\t', quote = TRUE
)
