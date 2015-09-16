

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
library(fpc)  # Density based clustering dbscan
library(gplots)  # Colorpanel
library(scatterplot3d)  # 3D plotting
library(monocle)
library(tsne)  # Non-linear ordination
library(pheatmap)
library(MASS)
library(cluster)
library(mclust)
library(flexmix)
library(lattice)
library(fpc)
library(amap)
library(RColorBrewer)
library(locfit)
library(Seurat)
```

Load Code
===
class:small-code


```r
# Source code
source(file.path("src","Modules.R")) # Helper functions
```

Briefly Single-cell Sequencing
===

#TODO
# UMI / RMT
# Sequence breakdown

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
data.set = read.delim( file.path("data","GSE29087_L139_updated.txt") )
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
rownames( data.set )
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
# ( how many genes express )
genes.per.cell <- apply( data.set, 2, function(x) sum( x>0 ))
```

How Many Expressed Genes (Complexity)?
===
class:small-code


```r
cell.outlier = plot.cell.complexity( genes.per.cell )
```

![plot of chunk unnamed-chunk-9](single_cell_analysis-figure/unnamed-chunk-9-1.png) 

Filter Cells: Removing the Outlier Cells
===
class:small-code

- Cells that are unusually simple (or no expression)
- Cells that are unusually complex


```r
cell.outlier
```

```
 ES_F06 MEF_D12 
     46      92 
```

```r
# Remove outlier cells
data.set = data.set[,-1*cell.outlier]
ncol( data.set)
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

Genes Have Different Distributions
===
class:small-code


```r
plot.quantiles( data.set ) 
```

![plot of chunk unnamed-chunk-11](single_cell_analysis-figure/unnamed-chunk-11-1.png) 

Filter Genes: Using Prevalence
===
class:small-code


```r
# Remove low expressing genes
data.cleaned <- func_filter_by_occurence( data.set, 10, 10 )
```

Normalization in scData
===

- Lack of publications / annecdotal .
- Log( Count / cell sum * median magnitude ) + 1 .
- Median magnitude = Median of cell medians .

Normalizing for Cell Sequencing Depth
===
class:small-code


```r
data.cleaned.norm <- func_cpx( data.cleaned )
```
- Some perform TMM normalization afterwards.
  - See EdgeR package.

Sequencing Saturation
====

- The correct depth of sequencing will depend on the cell and the question.
- Can view saturation levels.

Plotting Sequencing Saturation of a Cell
===
class:small-code


```r
func_plot_saturation_curve( data.cleaned[,1], 1000 )
```

![plot of chunk unnamed-chunk-14](single_cell_analysis-figure/unnamed-chunk-14-1.png) 

Quality Control in scData
===

#TODO Izzary paper

Loading Data into Seurat
===
class:small-code


```r
nbt=read.into.seurat(file.path("data","HiSeq301_RSEM_linear_values.txt"),sep="\t",header=TRUE,row.names=1)
nbt=setup(nbt,project="NBT",min.cells = 3,names.field = 2,names.delim = "_",min.genes = 1000,is.expr=1)
```

Quality Control in scData
===
class:small-code

- Check the identity of the cells!!!

```r
vlnPlot( nbt, c("DPPA4") )
```

Viewing Specific Genes in Data
===
class:small-code

- Check the identity of the cells!!!

![plot of chunk unnamed-chunk-17](single_cell_analysis-figure/unnamed-chunk-17-1.png) 

Viewing Genes vs Genes
===
class:small-code


```r
cellPlot(nbt,nbt@cell.names[1],nbt@cell.names[2],do.ident = FALSE)
```

Viewing Genes vs Genes
===
class:small-code

![plot of chunk unnamed-chunk-19](single_cell_analysis-figure/unnamed-chunk-19-1.png) 

Viewing Cells vs Cells
===
class:small-code


```r
cellPlot(nbt,nbt@cell.names[3],nbt@cell.names[4],do.ident = FALSE)
```

Viewing Cells vs Cells
===
class:small-code

![plot of chunk unnamed-chunk-21](single_cell_analysis-figure/unnamed-chunk-21-1.png) 

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
nbt=prep.pca.seurat( y.cutoff = 2,x.low.cutoff = 2 )
pca.plot(nbt,1,2,pt.size = 3)
```

PCA using Seurat
===
class:small-code

![plot of chunk unnamed-chunk-23](single_cell_analysis-figure/unnamed-chunk-23-1.png) 

Identifying Genes Contributing to Components
===
class:small-code


```r
print.pca(nbt,1)
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

Identifying Genes Contributing to Components
===
class:small-code


```r
viz.pca(nbt,1:2)
```

![plot of chunk unnamed-chunk-25](single_cell_analysis-figure/unnamed-chunk-25-1.png) 

tSNE: What and Why?
===

tSNE using Seurat
===
class:small-code


```r
nbt=run_tsne(nbt,dims.use = 1:11,max_iter=2000)
tsne.plot(nbt,pt.size = 3)
```

tSNE using Seurat
===
class:small-code

![plot of chunk unnamed-chunk-27](single_cell_analysis-figure/unnamed-chunk-27-1.png) 

tSNE: PCA & tSNE side by side
===
class:small-code

![plot of chunk unnamed-chunk-28](single_cell_analysis-figure/unnamed-chunk-28-1.png) 

---

![plot of chunk unnamed-chunk-29](single_cell_analysis-figure/unnamed-chunk-29-1.png) 

QC the Clusters!
===
class:small-code


```r
pca.plot(nbt,1,2,pt.size = 3,group.by="nGene")
tsne.plot(nbt,pt.size = 3,group.by="nGene")
```

QC the Clusters!
===
class:small-code

![plot of chunk unnamed-chunk-31](single_cell_analysis-figure/unnamed-chunk-31-1.png) 

---

![plot of chunk unnamed-chunk-32](single_cell_analysis-figure/unnamed-chunk-32-1.png) 

Defining Clusters from PCA or TSNE
===

- Density based clustering.

```r
nbt=DBclust_dimension(nbt,1,2,reduction.use="tsne", G.use=8,set.ident=TRUE)
nbt=buildClusterTree(nbt,do.reorder=TRUE,reorder.numeric=TRUE,pcs.use=1:11,do.plot=FALSE)
tsne.plot(nbt,do.label=TRUE,label.pt.size=0.5)
```

Defining Clusters from PCA or TSNE
===

![plot of chunk unnamed-chunk-34](single_cell_analysis-figure/unnamed-chunk-34-1.png) 

Seurat: What Did We Miss?
===

- Data imputation
- Spatial Inference

Differential Expression in Defined Groups
===

- Biomarker discovery.
- What are current options?
 - ROC / T-test / LRT / Tobit-censored LRT tests ( Seurat ).
 - ( Monocle ).
 - SCDE ( Single-Cell Differential Expression ).

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
## Get groupings
#data.groups <- rep( NA, ncol( data.cleaned ) )
#data.groups[ grep( "MEF", names( data.cleaned )) ] <- "MEF"
#data.groups[ grep( "ES", names( data.cleaned )) ] <- "ES"
#data.groups <- factor( data.groups, levels = c("ES","MEF") )
```

SCDE: in Code
===
class:small-code


```r
#library(scde)

## Calculate error models
#o.ifm <- scde.error.models( as.matrix( data.cleaned ), groups = data.groups, n.cores=3, threshold.segmentation=TRUE, save.crossfit.plot=FALSE, save.model.plots=FALSE, verbose=1 )

## Filter out cell (QC)
#o.ifm <- o.ifm[ o.ifm$corr.a > 0, ]
```

SCDE: in Code
===
class:small-code


```r
## Set up the Prior (starting value)
#o.prior <- scde.expression.prior(models=o.ifm,counts=as.matrix( data.cleaned ), length.out=400,show.plot=FALSE)

## Perform T-test like analysis
#ediff <- scde.expression.difference(o.ifm,as.matrix(data.cleaned), o.prior,groups=data.groups,n.randomizations=100, n.cores=1,verbose=1)
#write.table(ediff[order(abs(ediff$Z),decreasing=T),], file="scde_results.txt",row.names=T,col.names=T, sep="\t",quote=F)
```

Visualize Differentially Expressed Genes
===
class:small-code


```r
## Read in results
#scde.results.de <- read.delim( "scde_results.txt", row.names=1 )
#head( scde.results.de )
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

Visualize Differentially Expressed Genes
===
class:small-code

Let's plot the top 100 DE genes.
 

```r
## Get top genes
#top.de.genes <- rownames( scde.results.de )[ 1:100 ]
 
## Visualize in heatmap
#heatmap( log2( data.cleaned[ top.de.genes,] + 1 ), vctr_grouping=data.groups )
```
 
Visualize Differentially Expressed Genes
===
class:small-code


```r
#top.de.genes <- rownames( scde.results.de )[ 1:100 ]
#heatmap( log2( data.cleaned[ top.de.genes,] + 1 ), vctr_grouping = data.groups )
```

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
source(file.path("src","RaceID_class.R"))
# Load tutorial data
race.in <-read.csv(file.path("data","transcript_counts_intestine.xls"), sep="\t", header=TRUE)
rownames( race.in ) <-race.in$GENEID
race.in <- race.in[grep("ERCC",rownames(race.in),invert=TRUE),-1]
race.data <- SCseq(race.in)
race.data <- filterdata( race.data, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000 )
race.data <- clustexp( race.data, metric="pearson", cln=0, do.gap=TRUE, clustnr=20, B.gap=50, SE.method="Tibs2001SEmax", SE.factor=.25, bootnr=50, rseed=17000 )
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
plotgap(race.data)
```

![plot of chunk unnamed-chunk-42](single_cell_analysis-figure/unnamed-chunk-42-1.png) 

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
race.data <- findoutliers( race.data, outminc=5, outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75)
race.data <- comptsne( race.data,rseed=15555 )
```

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
plottsne( race.data, final=FALSE )
```

![plot of chunk unnamed-chunk-44](single_cell_analysis-figure/unnamed-chunk-44-1.png) 

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
plottsne( race.data, final=TRUE )
```

![plot of chunk unnamed-chunk-45](single_cell_analysis-figure/unnamed-chunk-45-1.png) 

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
target.genes <- c("Apoa1__chr9", "Apoa1bp__chr3", "Apoa2__chr1", "Apoa4__chr9", "Apoa5__chr9")
plotexptsne( race.data, target.genes )
```

![plot of chunk unnamed-chunk-46](single_cell_analysis-figure/unnamed-chunk-46-1.png) 

RaceID: Detecting Rare Cell Populations
===
class:small-code


```r
plotsymbolstsne( race.data, type=sub("\\_\\d+$","",names(race.data@ndata)))
```

![plot of chunk unnamed-chunk-47](single_cell_analysis-figure/unnamed-chunk-47-1.png) 

RaceIDs: Leftovers
===

- Alternative methods to evaluating clusters.
- Differential Expression between clusters.

Time-series Analysis: Monocle
===

TODO: the methodology briefly

What to Expect
===

TODO:Can we tie this into SCDE?
- Read in and prepare data.
- Filter and quality control.
- Select genes of interest.
 - Literature or differential expression analysis.
- Order cells by expression.
- Run pseudotemporal analysis.

Moncole's Assumptions
===

- Monocle's Assumptions.
 - Genes not splice variants.
 - Assumes a log normal distribution.
 - Does NOT normalize (library size, depth, technical batch).
 - Do NOT give it raw counts.

Read In and Format
===
class:small-code

Needed files -
- Expression file ( Genes (row) x Cells (col) ) .
- Cell Phenotype Metadata ( Cells (row) x Metadatum (col) ) .
- Gene Metadata ( Genes (row) x Medatatum (col) ) .


```r
# Do not run (For later)
# monocle.data <- make_cell_data_set( expression_file="monocle_exprs.txt", cell_phenotype_file="monocle_cell_meta.txt", gene_metadata_file="monocle_gene_meta.txt" )

# Get data for today
monocle.data <- get_monocle_presentation_data()
```

Filter Genes
===
class:small-code


```r
# Require a minimun of 0.1 expression
monocle.data <- detectGenes( monocle.data, min_expr=0.1 )

# Require atleast 50 cells to have the minimum 0.1 expression
# Get name of genes pass these filters
monocle.expr.genes <- row.names(subset(fData( monocle.data), num_cells_expressed >= 50 ))
```

Other QC can be performed.
 - Depth, accurate capture of 1 cell, ...
 - Can be added to phenotype or feature metadata files.
 
Confirm the Log-Normal Assumption
===
class:small-code


```r
plot_log_normal_monocle( monocle.data )
```

![plot of chunk unnamed-chunk-50](single_cell_analysis-figure/unnamed-chunk-50-1.png) 

Ordering by Expression: Study View
===
class:small-code


```r
# Marker genes of biological interest
marker.genes <- get_monocle_presentation_marker_genes()
# Select from those marker genes those important to the study
ordering.genes <- select_ordering_genes( monocle.data, monocle.expr.genes, marker.genes, "expression~Media", 0.01 )

# Order the cells by expression
monocle.data <- order_cells_wrapper( monocle.data, ordering.genes, use_irlba=FALSE, num_paths=2, reverse=TRUE)
# Plot all cells in study with ordering
plot_spanning_tree( monocle.data )
# Write data to file
# TODO
```

Ordering by Expression: Study View
===
class:small-code


```
<simpleError in lm.fit(X.vlm, y = z.vlm, ...): NA/NaN/Inf in 'y'>
```

![plot of chunk unnamed-chunk-52](single_cell_analysis-figure/unnamed-chunk-52-1.png) 

Ordering by Expression: Gene View
===
class:small-code


```r
monocle.data.diff.states <- monocle.data[ monocle.expr.genes, pData( monocle.data)$State != 3]
# TODO subset is not generic, needs gene_short_name attr
subset.for.plot <- subset_to_genes( monocle.data.diff.states, c("CDK1","MEF2C", "MYH3") )
plot_genes_in_pseudotime( subset.for.plot, color_by="Hours" )
```

Ordering by Expression: Gene View
===
class:small-code

![plot of chunk unnamed-chunk-54](single_cell_analysis-figure/unnamed-chunk-54-1.png) 

Genes Which Follow an Assumed Temporal Pattern
===
class:small-code


```r
# Get genes of interest
subset.pseudo <- subset_to_genes( monocle.data, c("MYH3","MEF2C","CCNB2","TNNT1"))
subset.pseudo <- subset.pseudo[,pData(subset.pseudo)$State != 3 ]
subset.pseudo.diff <- differentialGeneTest( subset.pseudo, fullModelFormulaStr="expression~sm.ns(Pseudotime)")
plot_genes_in_pseudotime( subset.pseudo, color_by="Hours")
```

Genes Which Follow an Assumed Temporal Pattern
===
class:small-code

![plot of chunk unnamed-chunk-56](single_cell_analysis-figure/unnamed-chunk-56-1.png) 

Monocle Leftovers
===

- Other methodology.
 - Simple Differential Expression.
     - Finding genes that distiguish cell state.
 - Multifactorial Differential Expression.
     - Account for confounders.
 - Clustering genes by pseudotime.
 
Diffusion Maps
===

#TODO Add theory about diffusion maps.

What Did We Find ?
===

- Unsupervised discovery of substructure.
- Inference on known structure.
- How do we connect this to biology?

Where Do We Go?
===

- Gene Set Enrichment.
  - DAVID (online or R library RDAVIDWebService).
  - GSEA (online or many libraries).
    - wilcoxGST from the limma library.
    - GSEABase.
- GenePattern workshop.

Summary: of the Data
===

- We are still understanding scData and how to apply it.
  - Not normal.
  - Zero-inflated.
  - Very noisey.
  - Vary in library complexity.
- Keeping these characteristics in analysis assumptions.

Summary: of Today
===

- Created expectations for scData.
- Performed QC.
- Plotted genes and cells.
- Detected novel structure.
- Applied a statistical inference method.
- Performed pseudotemporal analysis.

Thank You
===

- Aviv Regev
- Asma bankapur
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

Notes: to Make a PDF
===
class:small-code

- Create a pdf file before you plot ( can plot multiple plots ).
- Close the plotting.


```r
#pdf( "data/my_file.pdf", useDingbats = FALSE ) # Start pdf
#plot( 1:10, log(1:10 ) ) # plot in to the pdf file
#plot( seq(0,.9,.1), sin(0:9) ) # another plot for the pdf file
#dev.off() # Close pdf file ( very important )
```

Mixture Modeling to Select Cells
===

- Uses an EM algorithm optimizing the number of gaussian distributions.
- Works directly off of the ordination.
- Only as good as the ordination.
- Cell ordination shows the strongest signals.
- Really easy!

---

![mclust_description](images/mclust.pdf)

mclust: Mixture Modeling
===
class:small-code


```r
#library(mclust) # Load library
## Start with our first two dimensions
#mclust.results = Mclust(results.pca$rotation[,c(1:2)])

## Get classification groups
#mclust.groups = mclust.results$classification

## Plot
#plot( mclust.results, what=c("classification") )
```

mclust: Mixture Modeling
===
class:small-code



Supervised vs Unsupervised
===
class:small-code



---


