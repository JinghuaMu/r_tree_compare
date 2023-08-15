## Description
The purpose of this repository is to conduct a horizontal comparison of the effects of single model and mixed model phylogenetic trees of iqtree, using realistic datasets and statistical methods for analysis.

## Composition
- `single_gene_compare.html`: The main script for conducting a horizontal comparison of the models.
- `iqtree_info.r`: Extracts the runtime results and parameter information from `.iqtree` files. Load this script by using `source("iqtree_info.r")`. Extract the information list using `summarise_iqtree("iqtree_file")` and reference it using the `$` operator. Use `short_description()` on the resulting list to generate quick information.
