# taxPPro 1.0.2

* When sources of the same confidence have a value for a taxon,
  `getSetWithIDs()` and `getSetWithoutIDs()` keep the source that comes
  first in the `Priority` column of bugphyzz's `attribute_sources.tsv`
  (bugphyzz >= 1.9.2), through the new `sourcePriority()`, in place of the
  alphabetically last source name (waldronlab/bugphyzzExports#66). Numeric
  attributes are not affected.

# taxPPro 1.0.0

* This first version was used for the v1.0.0 relase of bugphyzzExports.
