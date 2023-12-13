# taxPPro - Taxonomic and Phylogenetic Propagation

A provitional package to hold functions for ASR and inheritance
for bugphyzz


The package can be installed with:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install('sdgamboa/taxPPro')
```


Types of data (spreadsheets)

**Binary data with only TRUE values -- binary-T**
+ Merge Attribute and Atribute_value.
+ Don't propagate.
 
**Binary data with TRUE and FALSE values -- binary-TF**
+ Merge Attribute and Attribute_value.
+ Complete binary data with 0s.
+ Assign 0.5 to unknown tips.

**Multistate-Union with binary data with only TRUE values -- multistate-union-T**
+ Merge Attribute and Attribute_value.
+ Don't propagate.

**Multistate-Union with binary data with both TRUE and FALSE values -- multistate-union-TF**
+ Separate.
+ Merge Attribute and Attribute_value.
+ Complete binary data with 0s.
+ Assign 0.5 to unknown tips.

**Multistate-intersection -- multistate-intersection**
+ Completey data with 0s.
+ Assign 1 / n attributes to unknown tips.
+ Do not export inputed zeros (when data is completed, after binning, and after inheritance).

**Numeric -- numeric**
+ Convert to range.

**Range -- range**
+ Convert to categorical.
