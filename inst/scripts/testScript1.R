# library(bugphyzz)
# library(purrr)
# library(dplyr)
#
# phys <- physiologies(keyword = 'all', remove_false = FALSE, full_source = TRUE)
#
# binaries <- keep(phys, ~ length(unique(.x$Attribute)) == 1 && !'Attribute_value_max' %in% colnames(.x))
# ranges <- keep(phys, ~ length(unique(.x$Attribute)) == 1 && 'Attribute_value_max' %in% colnames(.x))
# categorical <- phys[!names(phys) %in% c(names(binaries), names(ranges))]
