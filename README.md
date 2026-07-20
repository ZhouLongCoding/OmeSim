# OmeSim

OmeSim is a Java-based simulator for generating molecular omics data together
with multiple phenotypic traits. In its first release, OmeSim simulates
whole-transcriptome expression data and whole-phenome trait data from
user-provided genotype, gene annotation, pathway, and parameter files.

OmeSim is designed for evaluating statistical models and computational tools
that integrate omics data to study the genetic basis of complex traits. It can
simulate additive and nonlinear genetic architectures, including epistasis,
compensatory, heterogeneous, and compound models, and it outputs both simulated
data and "gold-standard" correlation, association, and causality information.

Authors: Caifeng Li*, Zhou Long*, and Qingrun Zhang.

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Run OmeSim](#run-omesim)
- [Input files](#input-files)
- [Parameter file](#parameter-file)
- [Output files](#output-files)
- [Causality query](#causality-query)
- [Contact](#contact)
- [License](#license)

## Requirements

### Java

OmeSim is distributed as an executable JAR file. A working Java Runtime
Environment or Java Development Kit is required.

Check that Java is available:

```bash
java -version
```

If Java is not found, install Java and make sure the `java` command is available
in your shell `PATH`.

### R and R packages

R is required for OmeSim visualization outputs. Check that `Rscript` is
available:

```bash
which Rscript
```

If the default `Rscript` path does not work, set the correct path in
`parameter.txt` using `Rscript_Binary_Path`.

Install the required R packages:

```r
install.packages(c(
  "pheatmap",
  "tidyverse",
  "ggplotify",
  "heatmaply",
  "igraph",
  "ggraph"
))
```

## Installation

Download `OmeSim.tar.gz`, then decompress it:

```bash
tar -xvzf OmeSim.tar.gz
```

The package contains:

- `OmeSim.jar`: executable Java program
- `parameter.txt`: template parameter file
- `input_files`: example genotype, gene, and pathway data
- `visualization_scripts`: R scripts used to generate figures

## Run OmeSim

1. Decompress the package.
2. Edit `parameter.txt` so that all file paths match your local environment.
3. Run the simulation:

```bash
java -Xmx4g -jar OmeSim.jar Simulate -input parameter.txt
```

4. Check the folder specified by `Output_Folder`. In the examples below, this
   folder is named `results`.
5. If sample outputs are provided, compare the generated files with the sample
   outputs to confirm that the installation and configuration are working.

The `-Xmx4g` option gives Java up to 4 GB of memory. Increase this value if you
simulate large datasets and have enough memory available, for example
`-Xmx16g`.

## Input files

The `Simulate` function requires four main input files.

### 1. Genotype file

The genotype file provides genome-wide genetic variants. Supported formats are:

- `CSV`
- `PLINK`
- `VCF`

CSV is recommended for simple and efficient genotype storage. Genotypes should
be coded as:

- `0`: homozygous reference allele
- `1`: heterozygous allele
- `2`: homozygous alternative allele

Example CSV:

```csv
#Chr,Loc,Subject_1,Subject_2,Subject_3,Subject_4
1,34506,0,1,0,2
```

### 2. Gene file

The gene file defines genomic coordinates for genes. Columns are separated by
commas.

Example:

```csv
#Gene_ID,Chr,Start,End
ENSG00000000457.14,1,169849631,169894267
```

### 3. Pathway file

The pathway file defines pathway membership. Each line should include a pathway
ID, followed by a tab, followed by a comma-separated list of genes in that
pathway.

Example:

```text
Pathway_1<TAB>Gene_1,Gene_2,...,Gene_n
```

### 4. Parameter file

The parameter file controls input paths, output paths, genetic architecture,
simulation size, visualization, and model settings. A template `parameter.txt`
is included in the OmeSim package.

## Parameter file

At minimum, update the required paths in `parameter.txt`:

```ini
Genotype_File=/path/to/genotype_file.csv
Genotype_File_Format=CSV
Genes_File=/path/to/gene_file.csv
Pathway_File=/path/to/pathway_file.txt
Output_Folder=/path/to/results
Arch_Detailed_File=/path/to/causal/causal_terms.txt
Rscript_Binary_Path=/path/to/bin/Rscript
Visualization_Code_Folder=/path/to/visualization_scripts
```

Optional parameters can be added to the same file. If a parameter is not listed,
OmeSim uses its default value.

Example:

```ini
Num_Traits=1000
Min_MAF=0.05
Max_MAF=0.5
Loc_Distance=0
Iteration_Rounds=10
Gene_Contributors=1.0,0.5,0.3,0.1
Traits_Contributors=1.0,0.8,0.5,0.3
Gene_Contri_Weights=1.0,0.5,4.0,2.0
Traits_Contri_Weights=1.0,4.0,2.0
```

### Modify parameters

To change simulation settings, edit `parameter.txt`, save the file, and rerun
OmeSim.

For example, to simulate fewer traits:

```ini
Num_Traits=100
```

To change MAF filters:

```ini
Min_MAF=0.01
Max_MAF=0.5
```

To point OmeSim to a different R installation:

```ini
Rscript_Binary_Path=/usr/local/bin/Rscript
```

After editing the parameter file, rerun:

```bash
java -Xmx4g -jar OmeSim.jar Simulate -input parameter.txt
```

### Common parameters

| Parameter | Description | Default |
| --- | --- | --- |
| `Genotype_File` | Absolute path to the genotype file | Required |
| `Genotype_File_Format` | Genotype format: `CSV`, `PLINK`, or `VCF` | Required |
| `Genes_File` | Absolute path to the gene coordinate file | Required |
| `Pathway_File` | Absolute path to the pathway file | Required |
| `Output_Folder` | Absolute path to the main output directory, such as `/path/to/results` | Required |
| `Arch_Detailed_File` | Absolute path to `causal_terms.txt`, such as `/path/to/causal/causal_terms.txt`; `causal_terms.txt.tmp` is written to the same `causal` folder | Required |
| `Num_Traits` | Number of traits to simulate | `1000` |
| `Min_MAF` | Minimum minor allele frequency filter | `0.05` |
| `Max_MAF` | Maximum minor allele frequency filter | `0.5` |
| `Loc_Distance` | Minimum distance between adjacent variants | `0` |
| `Iteration_Rounds` | Number of iterative calculation rounds | `10` |
| `Gene_Contributors` | Probabilities for cis, trans, expression, and trait contributors to genes | `1.0,0.5,0.3,0.1` |
| `Traits_Contributors` | Probabilities for genetic, expression, trait, and infinitesimal contributors to traits | `1.0,0.8,0.5,0.3` |
| `Gene_Contri_Weights` | Average weights for contributors to gene expressions | `1.0,0.5,4.0,2.0` |
| `Traits_Contri_Weights` | Average weights for contributors to traits | `1.0,4.0,2.0` |
| `Weights_Relative_Range` | Relative fluctuation range around average contributor weights | `0.5` |
| `Traits_Var_Comp_Inf_Min` | Minimum variance component for the trait infinitesimal term | `0.05` |
| `Traits_Var_Comp_Inf_Max` | Maximum variance component for the trait infinitesimal term | `0.15` |
| `Var_Comp_Noise_Min` | Minimum variance component for the noise term | `0.05` |
| `Var_Comp_Noise_Max` | Maximum variance component for the noise term | `0.15` |
| `Cis_Var_Flanking` | Flanking region for cis-variants, in base pairs | `500000` |
| `Min_Var_Gene` | Minimum number of variants required for a gene | `10` |
| `Cis_Variant_Numbers_Mean` | Average number of cis-genetic variants contributing to a gene expression | `5` |
| `Cis_Variant_Numbers_Range` | Range around the mean number of cis-genetic variants contributing to a gene expression | `10` |
| `Trans_Variant_Numbers_Mean` | Average number of trans-genetic variants contributing to a gene expression; also controls trait genetic components | `5` |
| `Trans_Variant_Numbers_Range` | Range around the mean number of trans-genetic variants contributing to a gene expression; also controls trait genetic components | `10` |
| `Num_Contributing_Genes_Mean` | Average number of genes contributing to a gene expression or trait | `4` |
| `Num_Contributing_Genes_Range` | Range around the mean number of genes contributing to a gene expression or trait | `4` |
| `Num_Contributing_Traits_Mean` | Average number of traits contributing to a gene expression or trait | `2` |
| `Num_Contributing_Traits_Range` | Range around the mean number of traits contributing to a gene expression or trait | `1` |
| `Num_Infinitesimal` | Number of genome-wide variants selected to form the infinitesimal genetic background term | `5000` |
| `Gene_Models` | Probabilities for additive, epistatic, compensatory, heterogeneous, and compound gene models | `0.5,0.2,0.2,0.05,0.05` |
| `Trait_Models` | Probabilities for additive, epistatic, compensatory, heterogeneous, and compound trait models | `0.5,0.05,0.05,0.1,0.3` |
| `Trait_Binary_Proportion` | Proportion of binary traits | `0.5` |
| `Binary_Mode_Proportion` | Proportions of binary traits simulated by liability and logistic models | `0.5,0.5` |
| `Rscript_Binary_Path` | Absolute path to `Rscript` | `/usr/local/bin/Rscript` |
| `Visualization_Code_Folder` | Absolute path to the R visualization scripts | Required |
| `Large_Figures_Needed` | Whether to generate large full-dataset figures | `false` |
| `Trait_Causality_Degree` | Maximum degree for trait-specific sub-causality graphs | `4` |
| `Edge_Num_Ratio` | Size-control parameter for sub-causality graphs | `300` |

## Output files

OmeSim generates top-level result files in the `results` folder, a `figures` folder inside `results`, and a causal-term file in the `causal` folder.

Example output files are available [here](https://drive.google.com/drive/folders/1FXNOPCOzBTssAMDAuaLMLmBWHSY4l6gl). Users can compare their own generated files with these reference outputs to verify that OmeSim has been installed and executed correctly.

### Top-level result files

The main output files are:

- `asso_exp_by_genotype.csv`
- `asso_trait_by_genotype.csv`
- `causality_graph_one_many.csv`
- `causality_graph_one_one.csv`
- `corr_exp_by_exp.csv`
- `corr_trait_by_trait.csv`
- `corr_trait_by_exp.csv`
- `expressions.csv`
- `traits.csv`

In `expressions.csv` and `traits.csv`, each column represents an individual and
each row represents a gene expression or trait.

The correlation and association files record gold-standard relationships before
adding random noise and infinitesimal genetic background terms.

The causality graph files record relationships among genes and traits. Use
`causality_graph_one_many.csv` with the `Causality` function.

### Figures folder

When R and the required R packages are correctly configured, OmeSim also creates
a `figures` folder inside the `results` folder.

The `figures` folder contains one folder for each simulated trait and one folder
for each pathway. For the default setting of `Num_Traits=1000`, this means 1000
trait folders plus all pathway folders.

Each trait folder contains data and figures for the trait-specific causality
graph from degree 1 to degree 4.

Each pathway folder contains pathway-specific data and figures. For example,
the `Pathway_1` folder may contain:

- `Pathway_1.corr_exp_by_exp.csv`
- `Pathway_1.corr_exp_by_exp.png`
- `Pathway_1.expressions.csv`
- `Pathway_1.expressions.png`
- `Pathway_1.network.csv`
- `Pathway_1.network.png`
- `Pathway_1.trait_by_exp.csv`
- `Pathway_1.trait_by_exp.png`

### Causal folder

In addition to the files in `results`, OmeSim writes the following file to the
`causal` folder. This folder is defined by the path provided to
`Arch_Detailed_File`.

- `causal_terms.txt`

## Causality query

The `Causality` function checks causal relationships among queried term triples
after a simulation has been completed.

It can identify three relationship types:

- `Causality`: `variants_gene_ID => exp_gene_ID => trait_ID`
- `Pleiotropy`: `variants_gene_ID => exp_gene_ID` and
  `variants_gene_ID => trait_ID`
- `Reverse Causality`: `variants_gene_ID => trait_ID => exp_gene_ID`

Prepare a query file in which each line contains three comma-separated terms:

```csv
#variants_gene_ID,exp_gene_ID,trait_ID
ENSG00000278677.2,ENSG00000149480.7,Quantitative_Trait_0819
ENSG00000278677.2,ENSG00000149480.7,Quantitative_Trait_0000
ENSG00000169418.10,ENSG00000054282.16,Quantitative_Trait_0079
```

Run:

```bash
java -Xmx4g -jar OmeSim.jar Causality \
  -causal_graph_file /path/to/results/causality_graph_one_many.csv \
  -triple_IDs /path/to/query_terms.txt
```

The output is written automatically to:

```text
/path/to/query_terms.txt.checked.txt
```

Example output:

```text
causality: ENSG00000278677.2 => ENSG00000149480.7 => Quantitative_Trait_0819
pleiotropy: ENSG00000278677.2 => ENSG00000149480.7; ENSG00000278677.2 => Quantitative_Trait_0000
reverse_causality: ENSG00000169418.10 => Quantitative_Trait_0079 => ENSG00000054282.16
```

## Troubleshooting

### `java: command not found`

Install Java and confirm that `java -version` works in your terminal.

### R figures are not generated

Check that:

- R is installed.
- `Rscript_Binary_Path` points to the correct `Rscript` executable.
- `Visualization_Code_Folder` points to the folder containing the OmeSim R
  visualization scripts.
- Required R packages are installed.

### Out-of-memory errors

Increase the Java memory limit if your machine has enough available memory:

```bash
java -Xmx16g -jar OmeSim.jar Simulate -input parameter.txt
```

## Contact

Qingrun Zhang: <qingrun.zhang@ucalgary.ca>

## License

Copyright (c) 2026 Caifeng Li*, Zhou Long*, and Qingrun Zhang

OmeSim is released under the MIT Open Source License.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
