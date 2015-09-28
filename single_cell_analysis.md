

Single-cell RNA-Seq Analysis
========================================================
author: Timothy Tickle and Brian Haas
css: single_cell_analysis.css
date: October 1, 2015

Before We Get Started
===

- Single-cell analysis is new.
  - Give you a feel for the data.
  - Give you some options to explore.
  - These techniques will grow as the field does.

Before We Get Started
===

- This is VERY hands on.
  - Much can be applied to other analyses.
  - Strengthen those R ninja skills!
  - If you need, cut and pasting is available (cut_and_paste.txt).
  - Complex R is simplified in wrapper functions.
  
---

![ninja corgis](images/ninja_corgi.jpeg)

What We Will Attempt to Cover
===

- Initial data exploration, QC, and filtering.
- Ways to plot genes and cells.
- Finding clusters of cells.
- Performing differential expression.
- Detecting rare cell population.
- Psuedotemporal time-series analysis.

RStudio: Getting to Know You
===

Let's take a moment.
- Pull data from online (Github).
  - https://github.com/TimothyTickle/single_cell_analysis
- You can view this presentation on-line.
  - http://rpubs.com/timothyltickle/single_cell_analysis
- Quick overview of RStudio.

Logistics
===
class:small-code


```r
# Load libraries
library(vioplot) # violin plots
library(Seurat) # Data Basics, QC, Spatial
#library(scde) # Differential Expressions
#library(monocle) # Pseudotemporal analysis
```

Load Code
===
class:small-code


```r
# Source code Helper functions
source(file.path("src", "Modules.R"))
```

Briefly Single-cell Sequencing
===

[Drop Seq Video Abstract](http://www.sciencedirect.com/science/article/pii/S0092867415005498)

![drop_seq](images/dropseq_2.png)

Genes Have Different Distributions
===
class:small-code

![genes_log_normal](images/scdistributions_log_normal.png)

Genes Have Different Distributions
===
class:small-code

![genes_bimodal](images/scdistributions_bimodal.png)

Genes Have Different Distributions
===
class:small-code

![genes_zeroinfl](images/scdistributions_zeroinfl.png)

Genes Have Different Distributions
===
class:small-code

![genes_noise](images/scdistributions_noise.png)

Genes Have Different Distributions
===
class:small-code

![professor corgi](images/scdistributions_all.png)

Underlying Biology
===

Our First Data Set
===

Islam S et al. __Characterization of the single-cell transcriptional landscape by highly multiplex RNA-seq__ . Genome Research 2011

- 92 Cells.
  - Embryonic Stem Cells (ES).
  - Embryonic Fibroblasts (MEF).

Data: Ready, Get Set, Load!
===
class:small-code


```r
# Load tab delimited file
data.set = read.delim(file.path("data", "GSE29087_L139_updated.txt"))
```

Always Look at Your Data
===
class:midcenter

- These are important steps for any scRNA-Seq data set.

![professor corgi](images/professor_corgi.jpg)

What are Our Genes?
===
class:small-code


```r
# Gene names (row names)
rownames(data.set)
```

What are Our Genes?
===
class:small-code


```
[1] "Tor1aip2"      "Pnkd"          "Smyd3"         "4921521F21Rik"
[5] "Gpbar1"        "1700016C15Rik"
```

What are Our Cells?
===
class:small-code


```r
# Column names
# Sample / Cell names
colnames( data.set )
```

What are Our Cells?
===
class:small-code


```
[1] "ES_A01" "ES_B01" "ES_C01" "ES_D01" "ES_E01" "ES_F01"
```

How Many Expressed Genes (Complexity)?
===
class:small-code


```r
# Plot genes per cell
# How many genes express
genes.per.cell <- apply( data.set, 2, function(x) sum( x>0 ))
```

How Many Expressed Genes (Complexity)?
===
class:small-code


```r
# Each cell is a dot
cell.outlier = plot.cell.complexity( genes.per.cell )
```

<img src="single_cell_analysis-figure/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" style="display: block; margin: auto;" />

Filter Cells: Removing the Outlier Cells
===
class:small-code

- Cells that are unusually simple (or no expression)
- Cells that are unusually complex


```r
# Which outliers?
cell.outlier
```

```
 ES_F06 MEF_D12 
     46      92 
```

```r
# Remove outlier cells
data.set = data.set[, -1 * cell.outlier]
ncol(data.set)
```

```
[1] 90
```

Filter Cells: Removing the Outlier Cells
===
class:small-code

- Outlier samples are not just measured by complexity
  - Percent Reads Mapping
  - Percent Mitochondrial Reads
  - Presence of marker genes
  - etc ...

Filter Genes: Using Prevalence
===
class:small-code

![professor corgi](images/prob_ave_expression.png)

Filter Genes: Using Prevalence
===
class:small-code


```r
# Remove low expressing genes
# Remove genes that do not have atleast 10 counts in 10 samples.
data.cleaned <- func_filter_by_occurence( data.set, 10, 10 )
# Remove genes that their top 10 expression are not on average atleast 5
data.cleaned <- filter_by_average_top_gene( data.set, 10, 5 )
```

Normalization in scData
===

- Lack of publications / annecdotal .
- $CPX = Log( Count / cell sum * median\_magnitude ) + 1$.
- Median_magnitude = A magnitude more than the median of cell medians .
  - Cell Medians = 6024, Median magnitude = 10000
- Can also use TPX
  - Use RSEM for TPM
- Same for 3', same level of complexity

Normalizing for Cell Sequencing Depth
===
class:small-code


```r
# Counts to CPX
data.cleaned.norm <- func_cpx(data.cleaned)
```

Sequencing Saturation
====

- The correct depth of sequencing will depend on the cell and the question.
- Can view saturation levels.

Plotting Sequencing Saturation of a Cell
===
class:small-code


```r
func_plot_saturation_curve(data.cleaned[, 1], 1000)
```

<img src="single_cell_analysis-figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" style="display: block; margin: auto;" />

Summary: of the Data
===

- We are still understanding scData and how to apply it.
  - Not normal.
  - Zero-inflated.
  - Very noisey.
  - Vary in library complexity.
- Keeping these characteristics in analysis assumptions.

Seurat
===


Loading Data into Seurat
===
class:small-code


```r
# Read data in to Seurat
nbt = read.into.seurat(file.path("data", "HiSeq301_RSEM_linear_values.txt"), 
    sep = "\t", header = TRUE, row.names = 1)

# Sets up object, filters
nbt = setup(nbt, project = "NBT", min.cells = 3, names.field = 2, names.delim = "_", 
    min.genes = 1000, is.expr = 1)
```

Quality Control in scData
===
class:small-code

- Check the identity of the cells!!!

```r
# Plot gene expression throughout cell groups
vlnPlot(nbt, c("DPPA4"))
```

Viewing Specific Genes in Data
===
class:small-code

- Check the identity of the cells!!!

<img src="single_cell_analysis-figure/unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" style="display: block; margin: auto;" />

Viewing Genes vs Genes
===
class:small-code


```r
# Plot a gene vs a gene
cellPlot(nbt, nbt@cell.names[1], nbt@cell.names[2], do.ident = FALSE)
```

Viewing Genes vs Genes
===
class:small-code

<img src="single_cell_analysis-figure/unnamed-chunk-18-1.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" style="display: block; margin: auto;" />

Viewing Cells vs Cells
===
class:small-code


```r
# Plot a cell vs a cell
cellPlot(nbt, nbt@cell.names[3], nbt@cell.names[4], do.ident = FALSE)
```

Viewing Cells vs Cells
===
class:small-code

<img src="single_cell_analysis-figure/unnamed-chunk-20-1.png" title="plot of chunk unnamed-chunk-20" alt="plot of chunk unnamed-chunk-20" style="display: block; margin: auto;" />

Getting Your Data Out of R
===



Dimensionality Reduction and Ordination
===

- Start with many measurements (high dimensional).
  - Want to reduce to a few features (lower-dimensional space).
- One way is to extract features based on capturing groups of variance.
- Another could be to preferentially select some of the current features.
  - We have already done this.
- We need this to plot the cells in 2D (or ordinate them).

PCA: in Quick Theory
===

- Eigenvectors of covariance matrix.
- Find orthogonal groups of variance.
- Given from most to least variance.
  - Components of variation.
  - Linear combinations explaining the variance.
  
---

![pca_describe](images/PCA_plot_02.gif)

PCA: in Practice
===

Things to be aware of.
- Data with different magnitudes will dominate.
  - Zero center and divided by SD.
  - (Standardized).
- Can be affected by outliers.

---

![outlier_corgi](images/outlier_corgi.jpg)

PCA using Seurat
===
class:small-code


```r
# Prep data for PCA
nbt = prep.pca.seurat(y.cutoff = 2, x.low.cutoff = 2)

# Plot PCA
pca.plot(nbt, 1, 2, pt.size = 3)
```

PCA using Seurat
===
class:small-code

<img src="single_cell_analysis-figure/unnamed-chunk-22-1.png" title="plot of chunk unnamed-chunk-22" alt="plot of chunk unnamed-chunk-22" style="display: block; margin: auto;" />

Genes Contributing to Components
===
class:small-code


```r
# List driving PC1
print.pca(nbt, 1)
```

```
[1] "PC1"
 [1] "LGALS1"    "TIMP1"     "KRT18"     "IFI30"     "ARHGDIB"  
 [6] "IFI27"     "UCA1"      "HIST1H2BK" "KRT15"     "LCN2"     
[11] "S100A9"    "KRT81"     "ALDH1A3"   "KLK5"      "CEACAM6"  
[1] ""
 [1] "SOX11"     "TUBB2B"    "DCX"       "GPM6A"     "CRMP1"    
 [6] "RTN1"      "NNAT"      "C1orf61"   "STMN2"     "FABP7"    
[11] "LOC150568" "41520"     "TMSB15A"   "PPP2R2B"   "GAP43"    
[16] "NREP"     
[1] ""
[1] ""
```

Genes Contributing to Components
===
class:small-code


```r
# View genes driving PC1
viz.pca(nbt, 1:2)
```

<img src="single_cell_analysis-figure/unnamed-chunk-24-1.png" title="plot of chunk unnamed-chunk-24" alt="plot of chunk unnamed-chunk-24" style="display: block; margin: auto;" />

tSNE: Nonlinear Dimensional Reduction
===

![tsne](images/cell_cycle_1.png)

tSNE: Nonlinear Dimensional Reduction
===

![tsne](images/cell_cycle.png)

tSNE: What and Why?
===

![tsne](images/tsne.png)

tSNE using Seurat
===
class:small-code


```r
# Run the t-SNE alogirthm
nbt = run_tsne(nbt, dims.use = 1:11, max_iter = 2000)

# Plot results
tsne.plot(nbt, pt.size = 3)
```

tSNE using Seurat
===
class:small-code

<img src="single_cell_analysis-figure/unnamed-chunk-26-1.png" title="plot of chunk unnamed-chunk-26" alt="plot of chunk unnamed-chunk-26" style="display: block; margin: auto;" />

tSNE: PCA & tSNE side by side
===
class:small-code

![plot of chunk unnamed-chunk-27](single_cell_analysis-figure/unnamed-chunk-27-1.png) 

---

![plot of chunk unnamed-chunk-28](single_cell_analysis-figure/unnamed-chunk-28-1.png) 

Quality Control in scData
===

![batch_effects](images/batch_effects_2.png)

QC the Clusters!
===
class:small-code


```r
# Paint PCA by complexity
pca.plot(nbt, 1, 2, pt.size = 3, group.by = "nGene")

# Paint the t-SNE by complexity
tsne.plot(nbt, pt.size = 3, group.by = "nGene")
```

QC the Clusters!
===
class:small-code

![plot of chunk unnamed-chunk-30](single_cell_analysis-figure/unnamed-chunk-30-1.png) 

---

![plot of chunk unnamed-chunk-31](single_cell_analysis-figure/unnamed-chunk-31-1.png) 

Defining Clusters from PCA or TSNE
===

- Density based clustering.

```r
# Perform clustering
nbt = DBclust_dimension(nbt, 1, 2, reduction.use = "tsne", G.use = 8, set.ident = TRUE)

# Identify groups
nbt = buildClusterTree(nbt, do.reorder = TRUE, reorder.numeric = TRUE, pcs.use = 1:11, 
    do.plot = FALSE)

# Plot
tsne.plot(nbt, do.label = TRUE, label.pt.size = 0.5)
```

Defining Clusters from PCA or TSNE
===

<img src="single_cell_analysis-figure/unnamed-chunk-33-1.png" title="plot of chunk unnamed-chunk-33" alt="plot of chunk unnamed-chunk-33" style="display: block; margin: auto;" />

Getting Your Data Out of R
===



SCDE: Single Cell Differential Expression
===

![scde_intro](images/scde_intro.png)

SCDE: in Quick Theory
===

For each group (ES or MEF).
- Genes are modeled to have two groups of counts.
  - Noisey area highly prone to dropout (Poisson distribution).
  - "Amplified" signal (Negative Binomial distribution).
- This makes the error model or how much one can trust counts.
- Pairwise within groups.

Differential Expression.
- Expected value * the probability of dropout in that cell for that expression level .

SCDE: in Code
===
class:small-code


```r
## Setting up cells groups

# data.groups <- rep(NA, ncol(data.cleaned))
# data.groups[ grep("MEF", names(data.cleaned)) ] <- "MEF"
# data.groups[ grep("ES", names(data.cleaned)) ] <- "ES"
# data.groups <- factor(data.groups, levels = c("ES","MEF"))
```

SCDE: in Code
===
class:small-code


```r
#library(scde)

## Calculate error models (Takes time)
# o.ifm <- scde.error.models( as.matrix( data.cleaned ), groups = data.groups, n.cores=3, threshold.segmentation=TRUE, save.crossfit.plot=FALSE, save.model.plots=FALSE, verbose=1 )

## Filter out cell (QC)
# o.ifm <- o.ifm[ o.ifm$corr.a > 0, ]
```

SCDE: in Code
===
class:small-code


```r
## Set up the Prior (starting value)
# o.prior <- scde.expression.prior(models=o.ifm, counts=as.matrix(data.cleaned),length.out=400, show.plot=FALSE)

## Perform T-test like analysis
# ediff <- scde.expression.difference(o.ifm, as.matrix(data.cleaned), o.prior, groups=data.groups, n.randomizations=100, n.cores=1,verbose=1)

# write.table(ediff[order(abs(ediff$Z), decreasing=T),], file="scde_results.txt", row.names=T, col.names=T, sep="\t", quote=F)
```

Visualize Differentially Expressed Genes
===

- mle = log2 fold change (estimate) .
- ub and lb = upper and lower bound on mle .
- ce = log2 fold change (conservative estimate) .
- Z = Z-score .
- cZ = Z-score corrected for multiple hypothesis testing .

---

![scde_output](images/scde_output.png)

RaceID: Detecting Rare Cell Populations
===



RaceID: Detecting Rare Cell Populations
===
class:small-code

- Large clusters are identifed with K-means clustering
- Within each cluster a ploynomial is fit to measure outlier variance
  - Negative binomial
- New outlier clusters are made from the outlier cells

K-means Clustering
===

- Pick the number of clusters you want (K)
- Randomly pick K points to be the groups to start with
- Do the following until results do not change
  - Assign all other points to the nearest group
  - Update the group's center to the mean of the group

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
# Load code
source(file.path("src", "RaceID_class.R"))

# Load tutorial data
race.in <- read.csv(file.path("data", "transcript_counts_intestine.xls"), sep = "\t", 
    header = TRUE)

# Set up data
rownames(race.in) <- race.in$GENEID
race.in <- race.in[grep("ERCC", rownames(race.in), invert = TRUE), -1]
race.data <- SCseq(race.in)
```

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
# Filter data
race.data <- filterdata(race.data, mintotal = 3000, minexpr = 5, minnumber = 1, 
    maxexpr = 500, downsample = FALSE, dsn = 1, rseed = 17000)

# Cluster data
race.data <- clustexp(race.data, metric = "pearson", cln = 0, do.gap = TRUE, 
    clustnr = 20, B.gap = 50, SE.method = "Tibs2001SEmax", SE.factor = 0.25, 
    bootnr = 50, rseed = 17000)
```

```
boot 1 
boot 2 
boot 3 
boot 4 
boot 5 
boot 6 
boot 7 
boot 8 
boot 9 
boot 10 
boot 11 
boot 12 
boot 13 
boot 14 
boot 15 
boot 16 
boot 17 
boot 18 
boot 19 
boot 20 
boot 21 
boot 22 
boot 23 
boot 24 
boot 25 
boot 26 
boot 27 
boot 28 
boot 29 
boot 30 
boot 31 
boot 32 
boot 33 
boot 34 
boot 35 
boot 36 
boot 37 
boot 38 
boot 39 
boot 40 
boot 41 
boot 42 
boot 43 
boot 44 
boot 45 
boot 46 
boot 47 
boot 48 
boot 49 
boot 50 
```

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
# Check clusters
plotgap(race.data)
```

<img src="single_cell_analysis-figure/unnamed-chunk-39-1.png" title="plot of chunk unnamed-chunk-39" alt="plot of chunk unnamed-chunk-39" style="display: block; margin: auto;" />

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
# Find outliers
race.data <- findoutliers(race.data, outminc = 5, outlg = 2, probthr = 0.001, 
    thr = 2^-(1:40), outdistquant = 0.75)

# Project with t-SNE
race.data <- comptsne(race.data, rseed = 15555)
```

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
# Plot t-SNE
plottsne(race.data, final = FALSE)
```

<img src="single_cell_analysis-figure/unnamed-chunk-41-1.png" title="plot of chunk unnamed-chunk-41" alt="plot of chunk unnamed-chunk-41" style="display: block; margin: auto;" />

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
plottsne(race.data, final = TRUE)
```

<img src="single_cell_analysis-figure/unnamed-chunk-42-1.png" title="plot of chunk unnamed-chunk-42" alt="plot of chunk unnamed-chunk-42" style="display: block; margin: auto;" />

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
# Interesting genes
target.genes <- c("Apoa1__chr9", "Apoa1bp__chr3", "Apoa2__chr1", "Apoa4__chr9", 
    "Apoa5__chr9")

# Plot expresssion of interesting genes
plotexptsne(race.data, target.genes)
```

<img src="single_cell_analysis-figure/unnamed-chunk-43-1.png" title="plot of chunk unnamed-chunk-43" alt="plot of chunk unnamed-chunk-43" style="display: block; margin: auto;" />

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
# Plot t-SNE labeling groups with symbols
plotsymbolstsne(race.data, type = sub("\\_\\d+$", "", names(race.data@ndata)))
```

<img src="single_cell_analysis-figure/unnamed-chunk-44-1.png" title="plot of chunk unnamed-chunk-44" alt="plot of chunk unnamed-chunk-44" style="display: block; margin: auto;" />

Getting Your Data Out of R
===



Where Do We Go?
===

- Gene Set Enrichment.
  - DAVID (online or R library RDAVIDWebService).
  - GSEA (online or many libraries).
    - wilcoxGST from the limma library.
    - GSEABase.
- GenePattern workshop.

What Did We Miss?
===

- Seurat
  - Data imputation and Spatial Inference.
- SCDE
  - Batch Effect Correction.
- RaceIDs
  - Alternative methods to evaluating clusters.
  - Differential Expression between clusters.

Summary: of Today
===

- Created expectations for scData.
- Performed QC.
- Plotted genes and cells.
- Detected novel structure.
- Applied a statistical inference method.

Thank You
===

- Aviv Regev
- Asma Bankapur
- Brian Haas
- Itay Tirosh
- Karthik Shekhar
- Today's Attendees

References
===

Please note this is a collection of many peoples ideas.
Included in the download is a references.txt to document sources, tutorials, software, and links to cute corgi pictures :-)

Questions?
===

![gradute corgi](images/graduate_corgi.jpg)

Monocle
===


Time-series Analysis: Monocle
===

Orders expression into Pseudotime
  - Quantitative measurement of progress through a biological process
![monocle](images/monocle.png)

Moncole's Assumptions
===

- Monocle's Assumptions.
 - Genes not splice variants.
 - Assumes a log normal distribution.
 - Does NOT normalize (library size, depth, technical batch).
 - Do NOT use raw counts.

Read In and Format
===
class:small-code

Needed files -
- Expression file ( Genes (row) x Cells (col) ) .
- Cell Phenotype Metadata ( Cells (row) x Metadatum (col) ) .
- Gene Metadata ( Genes (row) x Medatatum (col) ) .





















```
Error in getClass(Class, where = topenv(parent.frame())) : 
  "AnnotatedDataFrame" is not a defined class
```
