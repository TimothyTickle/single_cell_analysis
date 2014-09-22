

scAnalysis
========================================================
author: Timothy Tickle and Brian Haas
css: 2014_scAnalysis.css
date: September 23, 2014

Before we get started
===

- scAnalysis is new
  - Give you a feel for the data
  - Give you some options to explore
  - These techniques will grow as the field does

Before we get started
===

- This is VERY hands on
  - Much can be applied to other analyses
  - Strengthen those R ninja skill!
  - If you need, cut and pasting is available
  
---

![ninja corgis](images/ninja_corgi.jpeg)

Before we get started
===

- The plotting shown here will not be ggplots
  - ggplots is AWESOME! but we are just trying to keep this simple

What we will *Attempt* to Cover
===

- How to look at the data *
- Ways to plot samples (ordination) *
- Discovery substructure in data *
- Using gene sets for visualization  *
- An example of unique plots for scAnalysis
- Discussion on performing inference

RStudio: getting to know you
===

Let take a moment
- Pull data from online (Github)
- Quick overview of RStudio

Logistics
===
class:small-code


```r
# Load libraries
library(boot)  #For SCDE
library(caret)  #Near-zero filter
library(GMD)  #Cluster selection on dendrograms
library(gplots)  #Colorpanel
library(mclust)  #Selection of clusters
library(scatterplot3d)  #3D plotting
library(scde)  #Single cell differential expression
library(vegan)  #NMDS and distance metrics

# Source code
source("heatmap.3b.R")  #Custom HCL
source("cell_cycle_plot.R")  #Custom plot for cell cycle
source("Modules.R")  #Helper functions
```

Today's data set
===

Islam S et al. __Characterization of the single-cell transcriptional landscape by highly multiplex RNA-seq__ . Genome Research 2011

- 96 Samples
  - Embryonic Stem Cells (ES)
    - 2.5 k distinct genes expressed
    - Found more correlated than MEF
  - Embryonic Fibroblasts (MEF)
    - 5.4 k distinct genes expressed
- Transcripts per Million (tpm)

Data: Ready, Start, Load!
===
class:small-code


```r
# Load tab delimited file
data = read.delim( file.path("data","GSE29087_L139_expression_tab.txt"), row.names = 1 )

# For convenience splitting in to metadata and data
metadata <- data[ 1:6 ]
data <- data[ -1 * 1:6 ]

# Remove features without counts
zero.features <- which( apply( data, 1, sum ) == 0 )
data <- data[ -1 * zero.features, ]
metadata <- metadata[ -1 * zero.features, ]
length( zero.features )
```

```
[1] 8023
```

Data: Ready, Start, Load!
===
class:small-code


```r
# Get groupings
data.groups <- c( rep("ES",48),rep("MEF",48))

# Get colors to plot with
data.groups.colors <- func_factor_to_metadata_color( as.factor( data.groups) )$vctr_grouping_colors
```
TODO update

Always look at your data
===

![professor corgi](images/professor_corgi.jpg)

A quick look at the data set
===
class:small-code


```r
# Get dimensions
dim( data )
```

```
[1] 14913    96
```

```r
# See first couple column names
names( data )[1:10]
```

```
 [1] "A01" "B01" "C01" "D01" "E01" "F01" "G01" "H01" "A02" "B02"
```

```r
# Summary per column
summary(data)
```

A quick look at the data set
===
class:small-code


```r
dim( metadata )
```

```
[1] 14913     6
```

```r
names( metadata )
```

```
[1] "Chr"         "Pos"         "Strand"      "TrLen"       "MinExonHits"
[6] "MaxExonHits"
```

```r
summary(metadata)
```

Let's characterize the data
===

- Is the data normal?
- Sparsity / zero-inflation
- Overdispersed

Gene sparsity
===
class:small-code


```r
# Sum each gene
feature.sum <- apply( data, 1, sum )

# Sort sums
feature.sum.sorted <- sort( feature.sum )

# Plot with quartiles
plot( log( feature.sum.sorted ), main = "Total gene counts throughout samples", ylab = "Log total counts", xlab = "Genes in sorted order" )

# Add quartile lines
# Min 1, 1st quartile 3729, Median 7457, 3rd quartile 11180, Max 14913
abline( v = c(1, 3729,7457,11180,14913), col = c("red", "violet","cyan","blue","green"))
```

Gene sparsity
===
class:midcenter

![plot of chunk unnamed-chunk-9](scAnalysis-figure/unnamed-chunk-9.png) 

Normality?: at different sparsity levels
===
class:small-code


```r
# Get the sorted order for the genes
feature.sum.order <- order( feature.sum )

# Get the 10 most sparse genes
index_min <- feature.sum.order[ 1:10 ]

# 1 quartile sparsity genes
index_q1 <- feature.sum.order[ 3724:3723 ]

# Median sparsity genes
index_median <- feature.sum.order[ 7452:7461 ]

# 3rd quartile sparsity genes
index_q3 <- feature.sum.order[ 11175:11184 ]

# Get the least sparse genes
index_max <- feature.sum.order[ 14908:14913 ]
```

Zooming in to gene distributions
===
class:small-code


```r
# Start plot
plot( x=0,y=0,type="p", xlim=c(0,100), ylim=c(0,.2) )

# Plot from least to most sparse genes
for( i_min_plot in index_min ){ lines( density( as.matrix( data[ i_min_plot, ])), col="red", add=TRUE)}

for( i_q1_plot in index_q1 ){ lines( density( as.matrix( data[ i_q1_plot, ])), col="violet", add=TRUE)}

for( i_median_plot in index_median ){ lines( density( as.matrix( data[ i_median_plot, ])), col="cyan", add=TRUE)}

for( i_q3_plot in index_q3 ){ lines( density( as.matrix( data[ i_q3_plot, ])), col="blue", add=TRUE)}

for( i_max_plot in index_max ){ lines( density( as.matrix( data[ i_max_plot, ])), col="green", add=TRUE)}
```

Zooming in to gene distributions
===

Many zeros
![plot of chunk unnamed-chunk-12](scAnalysis-figure/unnamed-chunk-12.png) 

---

Zoomed in
![plot of chunk unnamed-chunk-13](scAnalysis-figure/unnamed-chunk-13.png) 

Zooming in to gene distributions
===

Log transform?
![plot of chunk unnamed-chunk-14](scAnalysis-figure/unnamed-chunk-14.png) 

---

Zoomed in 
![plot of chunk unnamed-chunk-15](scAnalysis-figure/unnamed-chunk-15.png) 

Sparsity
===
class:small-code


```r
# Check samples for read depth
sample.depth <- apply( data, 2, sum )
summary(sample.depth )
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  10100  183000  337000  556000  941000 2220000 
```

Sparsity
===
class:small-code


```r
barplot( sort( sample.depth ), main = "Sample Depth ", xlab = "Sample", ylab = "Depth")
```

Sparsity
===
class:midcenter

![plot of chunk unnamed-chunk-18](scAnalysis-figure/unnamed-chunk-18.png) 

Sparsity
===
class:small-code


```r
# Get Percent zero per gene
feature.percent.zero <- apply( data, 1, function(x){ length( which(x==0))/length(x)})

# Get mean per gene (excluding zeros)
feature.mean.no.zero <- apply( data, 1, function(x){ mean( x[ which(x!=0) ] ) })

#Plot with loess line
plot( feature.percent.zero, log( feature.mean.no.zero ), xlab = "log mean", ylab = "percent zero", main = "Feature Sparsity by Expression" )

lines(loess.smooth( feature.percent.zero, log(feature.mean.no.zero),degree=2,span=.2),col="cyan")
```

Sparsity
===
class:midcenter

![plot of chunk unnamed-chunk-20](scAnalysis-figure/unnamed-chunk-20.png) 

Overdispersion
===
class:small-code


```r
# Get SD without zeros
feature.sd.no.zero <- apply( data, 1, function(x){ sd( x[ which(x !=0 )])})

# Remove NAs and zeros SD
remove.no.zero <- which( is.na( feature.sd.no.zero ))
remove.no.zero <- c( remove.no.zero, which( feature.sd.no.zero == 0 ))

# Remove NAs and zeros mean
remove.no.zero <- c( remove.no.zero, which( is.na( feature.mean.no.zero )))
remove.no.zero <- unique( c( remove.no.zero, which( feature.mean.no.zero == 0 )))

# Reduce the lists
feature.sd.no.zero <- feature.sd.no.zero[-1*remove.no.zero]
feature.mean.no.zero <- feature.mean.no.zero[-1*remove.no.zero]
```

Overdispersion
===
class:midcenter

![plot of chunk unnamed-chunk-22](scAnalysis-figure/unnamed-chunk-22.png) 

Can QC Help?
===

- Removing very sparse features
- Imputing outliers

Removing Sparse Features
===
class:small-code


```r
# Get the percentile at 50% for each sample
sample.percentile <- apply( data, 2, function(x){ quantile(x[x !=0 ], .5)})
# Select noise features
feature.noise <- which(apply( data, 1, min_occurence_at_min_value, sample.percentile ) <= 10)
# Order the feauture by their total expression
feature.noise.by.expression <- order( apply( data[ feature.noise, ], 1, sum ), decreasing = TRUE)
```

Removing Sparse Features
===
class:small-code


```r
plot(density( as.matrix(data[ feature.noise[ feature.noise.by.expression[ 1 ]], ])))
```

Removing Sparse Features
===

![plot of chunk unnamed-chunk-25](scAnalysis-figure/unnamed-chunk-25.png) 

---

![plot of chunk unnamed-chunk-26](scAnalysis-figure/unnamed-chunk-26.png) 

Removing Sparse Features
===
class:small-code


```r
dim( data )
```

```
[1] 14913    96
```

```r
data <- data[ -1 * feature.noise, ]
metadata <- metadata[ -1 * feature.noise, ]
dim( data )
```

```
[1] 5654   96
```

Sample Read Depth: revisited
===
class:small-code


```r
# Recalculate sample depth
sample.depth <- apply( data, 2, sum )

# Make a color palette based on depth
depth.colors <- polychromatic_palette( length( sample.depth ) )

# Summarize sample depth
summary( sample.depth )
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  10000  177000  321000  522000  878000 2080000 
```

Sample Read Depth: revisited
===
class:small-code


```r
barplot( sort(sample.depth ) )
```

Sample Read Depth: revisited
===
class:midcenter

![plot of chunk unnamed-chunk-30](scAnalysis-figure/unnamed-chunk-30.png) 

Handling Challenges in the Data
===

- Can Transforms / Normalization Help?
  - Log transform seemed to work for less sparse
    - High expression
- We have a large range of counts per sample
  - In today's dataset this may be biology...
- If you want to remove them
  - Ordination, removing the dimension most correlated with depth
  - Inference, controlling for depth
- Rarifying?
- __Waste Not Want Not__ PLOS Comp. Bio. 2014

Dimensionality Reduction and Ordination
===

- Start with many measurements (high dimensional)
  - Want to reduce to a few features (lower-dimensional space)
- One way is to extract features based on capturing groups of variance
- Another could be to preferentially select some of the cuurent features
  - We have already done this
- We need this to plot the samples in 2D (or ordinate them)

PCA: in quick theory
===

- Eigenvectors of covariance matrix
- Find orthognal groups of variance
- Given from most to least variance
  - Components of variation
  - Linear combination of variance

PCA: in practice
===

Things to be aware of
- Data with different magnitudes will dominate
  - Zero center and divided by SD
  - (Standardized)
- Can be affected by outliers

PCA: in code
===
class:small-code


```r
# Row center and log
data.scaled <- t( scale( t( as.matrix( log( data + 1 ) ) ), center=TRUE, scale=TRUE ) )

# Remove constant rows
# data.scaled<-data.scaled[,-1*nearZeroVar(data.scaled)]

# Perfrom PCA
results.pca <- prcomp( data.scaled, retx = TRUE )
```

PCA: in code
===


```r
# Plot scree / elbow plot
plot( results.pca )
```

![plot of chunk unnamed-chunk-32](scAnalysis-figure/unnamed-chunk-32.png) 

PCA: in code
===
class:small-code


```r
# Get Percent variance
pca.var <- results.pca$sdev^2
pca.var <- pca.var/sum( pca.var )
pca.var <- round( pca.var, 2 )

# Scatter plot PCA
plot( results.pca$rotation[,1], results.pca$rotation[,2], pch=16, xlab=paste("PC1 (",pca.var[1],")"), ylab=paste("PC2 (",pca.var[2],")"), main="Standard PCA")
```

PCA: in code
===
class:midcenter

![plot of chunk unnamed-chunk-34](scAnalysis-figure/unnamed-chunk-34.png) 

PCA: by depth
===
class:small-code


```r
# Plot PCA
plot( results.pca$rotation[,1], results.pca$rotation[,2], pch=16, col=depth.colors, xlab="PC1", ylab="PC2", main="PCA by Sample Depth" )

# Add legend
legend( "topright", c("Max","Median","Min"), fill=c("yellow","black","purple") )
```

PCA: by depth
===
class:midcenter

![plot of chunk unnamed-chunk-36](scAnalysis-figure/unnamed-chunk-36.png) 

PCA: 3D with caution
===
class:midcenter

![plot of chunk unnamed-chunk-37](scAnalysis-figure/unnamed-chunk-37.png) 

PCA: 3D with caution
===
class:small-code


```r
# 3D plotting method for scatter plots
scatterplot3d( x=results.pca$rotation[,1], y=results.pca$rotation[,3], z=results.pca$rotation[,2], color = depth.colors, xlab="PC1", ylab="PC3", zlab="PC2", main="3D PCA using Components 1-3" )
```

Alternatives?
===

- PCA makes assumptions
  - Shows linear relationships
  - Many zeros / differing sparisty may be a component
- Nonmetric Multidimensional Scaling
- Weighted PCA
  - Different than Sparse PCA

NMDS: in quick theory
===

- Number of dimensions are selected
- Distance matrix is calculated
- Samples are positioned
- Distance betwen the rank order in the higher and lower dinmensions are calculated
  - Stress
- Positions are moved to be better
- This continues until stress is at a minimum

NMDS: in practice
===

![wizard corgi](images/wizard_corgi.jpg)

---

- The magic is in the metric
- By default you only get 2 dimensions

NMDS: in code (Bray-curtis)
===
class:small-code


```r
# Make proportional
nmds.data = sweep( data, 2, colSums( data ), "/")

# Perform NMDS with Bray Curtis distance metric
nmds.b.c.result = metaMDS( comm=t(nmds.data), distance="bray", k=2, autotransfer=FALSE, trymax=10)
```

NMDS: in code (Bray-curtis)
===
class:small-code


```
Run 0 stress 0.09202 
Run 1 stress 0.09317 
Run 2 stress 0.09418 
Run 3 stress 0.09531 
Run 4 stress 0.09742 
Run 5 stress 0.09531 
Run 6 stress 0.1533 
Run 7 stress 0.0963 
Run 8 stress 0.09821 
Run 9 stress 0.1026 
Run 10 stress 0.1911 
```

NMDS: in code (Bray-curtis)
===
class:small-code


```r
# Plot NMDS with Bray-curtis
plot( nmds.b.c.result$points[,1], nmds.b.c.result$points[,2], col=depth.colors, main="Ordination by Bray Curtis", xlab="NMDS 1", ylab="NMDS 2")

# Add legend
legend( "topleft", c("Max","Median","Min"), fill=c("yellow","black","purple") )
```

NMDS: in code (Bray-curtis)
===
class:midcenter

![plot of chunk unnamed-chunk-42](scAnalysis-figure/unnamed-chunk-42.png) 

NMDS: in code (Binomial)
===
class:small-code


```r
# Perform NMDS with binomial metric
nmds.b.result = metaMDS( comm=t(nmds.data), distance="binomial", k=2, autotransfer=FALSE, trymax=10)
```

```
Run 0 stress 0.0724 
Run 1 stress 0.09305 
Run 2 stress 0.09686 
Run 3 stress 0.09862 
Run 4 stress 0.08954 
Run 5 stress 0.09653 
Run 6 stress 0.094 
Run 7 stress 0.09774 
Run 8 stress 0.08516 
Run 9 stress 0.09001 
Run 10 stress 0.09692 
```

NMDS: in code (Binomial)
===
class:small-code


```r
# Plot and add legend
plot( nmds.b.result$points[,1], nmds.b.result$points[,2], col=depth.colors, main="Ordination by Binomial", xlab="NMDS 1", ylab="NMDS 2")

legend( "topleft", c("Max","Median","Min"), fill=c("yellow","black","purple") )
```

NMDS: in code (Binomial)
===
class:midcenter

![plot of chunk unnamed-chunk-45](scAnalysis-figure/unnamed-chunk-45.png) 

Next Steps in Ordination
===

Other options to explore
- Weighted PCA
- Distance metrics made for scAnalysis

Unsupervised Substructure Discovery
===

Often a goal of scProjects is to describe new structure to a group of cells:
- Heterogeniety of tumor populations
- Novel steps in development
- Robust / dynamic cellular signalling

Unsupervised Substructure Discovery
===

- The general approach
  - Select gene with most variance
  - Find cells that vary together (de novo cell groups)
  - Find the genes most likely to explain the cell groups
  - Further investigate the biology
    - Are these groupings important?

PCA + ANOVA: PCA
===

Overview
- Select major components of variation
- Group samples
- Find genes associated with those groups
- Group those genes

PCA + ANOVA: PCA
===
class:small-code

- Let's select major variance in the data set
- How many components do we select?


```r
#Plot scree plot
barplot( pca.var, xlab = "Components" )
number.components <- length( pca.var )
lines( 1:number.components, stick_breaking( number.components ), col="red" )
```

![plot of chunk unnamed-chunk-46](scAnalysis-figure/unnamed-chunk-46.png) 

PCA + ANOVA: PCA
===
class:small-code

 - Let's the to 1000 genes contributing to the variance
 

```r
# Find the order of the genes scores on component 1
ordered.pca.loadings <- order( results.pca$x[,1], decreasing=TRUE )

# Take the extreme 1000 genes
extreme.ends.pca <- c(ordered.pca.loadings[1:500], ordered.pca.loadings[5154:5654])
```

PCA + ANOVA: PCA
===
class:midcenter

![plot of chunk unnamed-chunk-48](scAnalysis-figure/unnamed-chunk-48.png) 


PCA + ANOVA: visualize
===
class:small-code

Let's see what the data looks like with a heatmap (HCL)


```r
# Subset the data to the extreme genes
data.scaled.subset <- data.scaled[ extreme.ends.pca, ]

# Visualize with HCL
heatmap( data.scaled.subset )
```

PCA + ANOVA: visualize
===
class:midcenter

![plot of chunk unnamed-chunk-51](scAnalysis-figure/unnamed-chunk-51.png) 



PCA + ANOVA: select sample groups
===
class:small-code


```r
# Create a distance matrix and cluster data
gmd.dist <- dist( t( data.scaled.subset ) )
hclust.aov <- hclust( gmd.dist )

# Evaluate clusters
gmd.clust <- css.hclust( gmd.dist, hclust.aov )
gmd.elbow.groups <- elbow.batch( gmd.clust, ev.thres=0.4, inc.thres=0.1 )
gmd.elbow.groups$k
```

```
[1] 3
```

```r
# Select groups
tree.groups <- as.factor( cutree( hclust.aov, k=gmd.elbow.groups$k ) )
```

PCA + ANOVA: select sample groups
===
class:small-code

These are the groups we have selected.


```r
heatmap( data.scaled.subset, vctr_grouping=tree.groups )
```

![plot of chunk unnamed-chunk-54](scAnalysis-figure/unnamed-chunk-54.png) 

PCA + ANOVA: PCA by Depth
===
class:midcenter

![plot of chunk unnamed-chunk-55](scAnalysis-figure/unnamed-chunk-55.png) 

PCA + ANOVA: Another view
===
class:small-code


```r
# Change groupings to colors
tree.group.colors = func_factor_to_metadata_color( tree.groups )$vctr_grouping_colors

# Plot with selectec groupings
plot( results.pca$rotation[,1], results.pca$rotation[,2], pch=16, xlab="PC1", ylab="PC2", main="PCA by Grouping (GMD)", col = tree.group.colors )
```

PCA + ANOVA: Another view
===
class:small-code

![plot of chunk unnamed-chunk-57](scAnalysis-figure/unnamed-chunk-57.png) 

Mixture modeling to select samples
===

- Uses an EM algorithm optimizing the number of guassian distributions
- Works directly off of the ordination
- Only as good as the ordination
- Sample ordination relys on the strongest signals
  - Global, less sensitive
- Really easy!

---

![mclust](images/mclust.pdf)

mclust: Mixture Modeling
===
class:small-code


```r
# Start with our first two dimensions
mclust.results = Mclust(results.pca$rotation[,c(1:2)])

# Get classification groups
mclust.groups = mclust.results$classification

# Plot
plot( mclust.results, what=c("classification") )
```

mclust: Mixture Modeling
===
class:midcenter

![plot of chunk unnamed-chunk-59](scAnalysis-figure/unnamed-chunk-59.png) 

PCA + ANOVA: select genes
===
class:small-code

Find genes signficantly associated with our sample groups


```r
# Create anova function to get p-values
do_aov <- function(x, tree.groups){
ret_aov = aov( x ~ tree.groups )
return( summary(ret_aov)[[1]][["Pr(>F)"]][1] ) }

# Call anova function on data frame
pca.list.p = apply( data.scaled.subset, 1, do_aov, tree.groups=tree.groups )

# Correct for multiple comparisons
pca.list.q = p.adjust( pca.list.p )

# Reduce to q-value cut off of 0.05
pca.list.names = names( pca.list.q[ pca.list.q <= 0.05 ] )
length( pca.list.names )
```

```
[1] 603
```

PCA + ANOVA: select Gene Groups
===
class:small-code


```r
# Reduce data frame to significant genes
data.scaled.subset = data.scaled.subset[ pca.list.names, ]

# Plot with heatmap (HCL)
heatmap( data.scaled.subset )
```

![plot of chunk unnamed-chunk-61](scAnalysis-figure/unnamed-chunk-61.png) 

PCA + ANOVA: select Gene Groups
===
class:small-code

We can select gene groups like we did sample groups


```r
# Get distance matrix and cluster data
gmd.dist.genes <- dist( data.scaled.subset )
hclust.aov.genes <- hclust( gmd.dist.genes )

# Evaluate clusters and select number of clusters
gmd.clust.genes <- css.hclust( gmd.dist.genes, hclust.aov.genes )
gmd.elbow.genes <- elbow.batch( gmd.clust.genes, ev.thres=0.4, inc.thres=0.1 )
gmd.elbow.genes$k
```

```
[1] 3
```

```r
# Break genes in to groups
tree.genes <- as.factor( cutree( hclust.aov.genes, k=gmd.elbow.genes$k ) )
```

Supervised Analysis: known gene groups
===

- Sometimes you may have very noisey data
- Our data set has strong signals
- Earlier it was mentioned one could average over gene groups
- Let's try an example!

Supervised Analysis: known gene groups
===
class:small-code


```r
# Read in data
gbm_data <- read.delim("data/GBM_data_matrix.txt", row.names=1)

# Reduce to genes of interest
gbm_CL <- read.delim("data/GBM_gene_CL.txt",stringsAsFactors=FALSE)
gbm_CL <- gbm_CL$genes
gbm_MES <- read.delim("data/GBM_gene_MES.txt",stringsAsFactors=FALSE)
gbm_MES <- gbm_MES$genes
gbm_NL<- read.delim("data/GBM_gene_NL.txt",stringsAsFactors=FALSE)
gbm_NL <- gbm_NL$genes
gbm_PN <- read.delim("data/GBM_gene_PN.txt",stringsAsFactors=FALSE)
gbm_PN <- gbm_PN$genes
```

Supervised Analysis: known gene groups
===
class:small-code


```r
# Reduce to genes of interest
gbm_all <- unique( c( gbm_CL, gbm_MES, gbm_NL, gbm_PN ) )
gbm_data <- gbm_data[ gbm_all, ]

# Visualize data initially
heatmap( as.matrix( gbm_data ) )
```

Supervised Analysis: known gene groups
===
class:midcenter

![plot of chunk unnamed-chunk-65](scAnalysis-figure/unnamed-chunk-65.png) 

Supervised Analysis: known gene groups
===
class:small-code


```r
# Make data fram to hold averages
gbm_averaged = as.data.frame( matrix(rep(NA,4*ncol( gbm_data )), ncol=ncol(gbm_data) ))
names( gbm_averaged ) <- names( gbm_data )
row.names( gbm_averaged ) <- c( "CL", "MES", "NL", "PN" )

# Add averages to data frame
gbm_averaged[1,] <- apply(gbm_data[ gbm_CL,],2,mean)
gbm_averaged[2,] <- apply(gbm_data[ gbm_MES,],2,mean)
gbm_averaged[3,] <- apply(gbm_data[ gbm_NL,],2,mean)
gbm_averaged[4,] <- apply(gbm_data[ gbm_PN,],2,mean)
```

Supervised Analysis: known gene groups
===
class:small-code


```r
# View after
heatmap( gbm_averaged )
```

![plot of chunk unnamed-chunk-67](scAnalysis-figure/unnamed-chunk-67.png) 

scOpportunities: cell cycle plot
===

Plot cell cycle with marker genes of interest

![cell cycle plot](images/cell_cycle.pdf)

Statistical Inference in scData
===

- Bayesian analysis
- Mixture Models
- Bimodal assumptions

SCDE: in quick theory
===

TODO

SCDE: in code
===
class:small-code


```r
# Setting up sampel groups
names( data.groups ) = data.groups
data.groups <- as.factor( data.groups )

# Calculate error models
o.ifm <- scde.error.models( as.matrix( data ), groups = data.groups, n.cores=1, threshold.segmentation=TRUE, save.crossfit.plot=FALSE, save.model.plots=FALSE, verbose=1 )
```

```
cross-fitting cells.
building individual error models.
WARNING: unable to find cross-fit models for  121  out of  1128  pairs. Using a subset.
WARNING: unable to find cross-fit models for  125  out of  1128  pairs. Using a subset.
```

```r
# Filter out cell (QC)
o.ifm <- o.ifm[ o.ifm$corr.a > 0, ]
```

SCDE: in code
===
class:small-code


```r
# Set up the Prior (starting value)
o.prior <- scde.expression.prior(models=o.ifm,counts=as.matrix(data), length.out=400,show.plot=FALSE)

# Perform T-test like analysis
ediff <- scde.expression.difference(o.ifm,as.matrix(data), o.prior,groups=as.factor(data.groups),n.randomizations=100, n.cores=1,verbose=1)
write.table(ediff[order(abs(ediff$Z),decreasing=T),], file="scde_results.txt",row.names=T,col.names=T, sep="\t",quote=F)
```

NMDS: in code (Reciprocal)
===
class:small-code


```r
# Make distance matrix
dist.recip <- reciprocal_weighting( x=data, groups=data.groups, n.cores=3)

# Perfroms NMDS
nmds.r.result = metaMDS( comm=dist.recip, k=2, autotransfer=FALSE )

#Plot
plot( nmds.r.result$points[,1], nmds.r.result$points[,2], col=data.groups.colors, main="Ordination by Reciprocal Weighting")
```

NMDS: in code (Reciprocal)
===
class:midcenter

![reciprocal](images/Reciprocal.pdf)

What did we find ?
===

- Unsupervised discovery of substructure
- Inference on known structure
- How do we connect this to biology?

Gene Set Enrichment Analysis: GSEA
===

- Not specifically related to scAnalysis
- Does help round out the story
- Covered in other section
- Many options
  - DAVID (online or R library RDAVIDWebService)
  - GSEA (online or many libraries)
    - wilcoxGST from the limma library
    - GSEABase
- GenePattern workshop

What did we find ???
===

- MEFs characterized by high expression of cytoskeletal proteins
- ESs characterized by highly enriched ribosomal proteins

TODO?

Summary: of the data
===

- We are still understanding scData and how to apply it
  - Not normal
  - Zero-inflated
  - Over-dispersed
  - Very noisey
- Keeping these characteristics in analysis assumptions

Summary: of today
===

- Created expectations for scData
- Applied 3 ordination techniques
- Tried 2 methods to detect substructure
- Applied known structure to a data set
- Applied 1 statistical inference method

Thank you
===

- Aviv Regev
- Alex Shalek
- Manik Kuchroo
- Rahul Satija
- Itay Tirosh

References
===

Please note this is a collection of many peoples ideas.
Included in the download is a references.txt to document sources.

Questions?
===

![gradute corgi](images/graduate_corgi.jpg)

Notes: to make a pdf
===

- Create a pdf file before you plot ( can plot multiple plots )
- Close the plotting


```r
pdf( "data/my_file.pdf", useDingbats = FALSE ) # Start pdf
plot( 1:10, log(1:10 ) ) # plot in to the pdf file
plot( seq(0,.9,.1), sin(0:9) ) # another plot for the pdf file
dev.off() # Close pdf file ( very important )
```

```
pdf 
  2 
```

