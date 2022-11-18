library(tibble)


df1 <- tribble(
    ~Taxon_name, ~Attribute, ~Score, ~Parent_taxon_name,
    'taxon1', 'aerobic', 1, 'Paredwnt',
    'taxon2', 'anaerobic', 1, 'same_parent',
    'taxon3', 'f_anaerobic', 1, 'same_parent'
)
