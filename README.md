## Description
The purpose of this repository is to conduct a horizontal comparison of the effects of single model and mixed model phylogenetic trees of iqtree, using realistic datasets and statistical methods for analysis.

## Composition
### Execuate script
- `single_gene_compare.rmd`: The main script for conducting a comparison of single and mixture models for one gene.
- `iqtree_info.r`: Extracts the runtime results and parameter information from `.iqtree` files. Load this script by using `source("iqtree_info.r")`. Extract the information list using `summarise_iqtree("iqtree_file")` and reference it using the `$` operator. Use `short_description()` on the resulting list to generate quick information.
- `dna_model.R`: Function dna_model takes a model string as input and calculates the rate matrix (Q) and other relevant parameters for DNA substitution models, returning a summary of the results. To use the function, pass a model string in .iqtree file as an argument. The output will include information about the model, the frequency type, the rate matrix (Q), and the base frequencies.
- `run_server.sh`: The script is a Bash script that processes multiple gene files using IQ-TREE software in parallel, saving the results in a single `.csv` file.
- `tree_comparison_combined.Rmd`: Improved version of the single gene alignment program script, accepts input from 'run server.sh' and automatically outputs the iqtree run result and gene tree alignment document in the specified directory.
### Analysis script
- `01_multi_compare.ipynb`: R script in jupyter format. Receive multiple csv results from `run_server.sh`.The program first merges multiple csv documents and organizes the data, adding new variables. Then perform information and visualize such as the category of the mixture model, the likelihood relative to the single model, the information criterion, the branch length, the running time (the csv file of the running time needs to be loaded separately).
- `02_try_modelling.ipynb`: R script in jupyter format. Accept merge statistics documents from '01_multi_compare.ipynb'. The program takes all the categorical variables and numerical variables obtained by single-class models as random variables, takes the class number of mixture model as predictor variables, fits linear models, and filters variables according to BIC criteria to determine the variables that are correlated with the class number of mixture model.
- `03_py_prediction.ipynb`: Python script in jupyter format. Accept merge statistics documents from '01_multi_compare.ipynb'. Functionally similar to the previous file, pycaret is used to determine the best machine learning model to predict the category of mixed models and to give variables with strong correlations.
- `04_inferred_tree.ipynb`: R script in jupyter format. Receiving species trees under multiple single-class models and mixed models, these species trees are compared by topological distance and qcf and posterior probability, and these species trees are compared with other species trees studied by scaly reptiles in the `treefile` folder.
- `05_subtitution_model.ipynb`: R script in jupyter format. Accept merge statistics documents from '01_multi_compare.ipynb'. Use a variety of statistics to analyze the characteristics and differences of different classes of alternative models in mixture models.
- `06_stabality.ipynb`: R script in jupyter format. Accept merge statistics documents from '01_multi_compare.ipynb'. Analyze whether the mixture model results in a more stable gene tree between different data sets.
### Other
- `coalescence_tree.sh`: Bash script. Merge all gene trees with specific names in the specified directory for subsequent use of astral to merge species trees.
- `Generate_UCE.R`: Generate fake UCE sequence for presentation.
- `time_record.py`: Python script. Merge the running time of gene tree analysis in all .iqtree documents in the specified directory.
### Data
All the generated Data is in the Data folder, and dragging and dropping the data into the root directory can directly execute scripts other than '01 multi compare.ipynb'.