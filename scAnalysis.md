

scAnalysis
========================================================
author: Timothy Tickle and Brian Haas
css: 2014_scAnalysis.css
date: September 23, 2014

Logistics
===


```r
library(caret)  #Near-zero filter
library(GMD)  #Cluster selection on dendrograms
library(gplots)  #Colorpanel
library(mclust)  # Selection of clusters
library(scde)  #Statistical inference (bayesian mixed model)
library(vegan)  #PCoA, distance metrics
source("heatmap.3b.R")
source("cell_cycle_plot.R")
source("Modules.R")
# RcppArmadillo, flexmix, multicore, Cairo
```

Today's Data Set
===

- Describe data set

How do we get it?
===


```r
# Load tab delimited file
data = read.delim( "data/GSE29087_L139_expression_tab.txt", row.names = 1 )

# For convenience splitting the data frame in to metadata and data
metadata = data[ 1:6 ]
data = data[ -1 * 1:6 ]

# Remove features without counts
zero.features <- which( apply( data, 1, sum ) == 0 )
length( zero.features )
```

```
[1] 8023
```

```r
data = data[ -1 * zero.features, ]
metadata = metadata[ -1 * zero.features, ]
```

Always look at your data
===

![Professor Corgi](images/professor_corgi.jpg)

A quick look at the data
===


```r
dim( data )
```

```
[1] 14913    96
```

```r
names( data )
```

```
 [1] "A01" "B01" "C01" "D01" "E01" "F01" "G01" "H01" "A02" "B02" "C02"
[12] "D02" "E02" "F02" "G02" "H02" "A03" "B03" "C03" "D03" "E03" "F03"
[23] "G03" "H03" "A04" "B04" "C04" "D04" "E04" "F04" "G04" "H04" "A05"
[34] "B05" "C05" "D05" "E05" "F05" "G05" "H05" "A06" "B06" "C06" "D06"
[45] "E06" "F06" "G06" "H06" "A07" "B07" "C07" "D07" "E07" "F07" "G07"
[56] "H07" "A08" "B08" "C08" "D08" "E08" "F08" "G08" "H08" "A09" "B09"
[67] "C09" "D09" "E09" "F09" "G09" "H09" "A10" "B10" "C10" "D10" "E10"
[78] "F10" "G10" "H10" "A11" "B11" "C11" "D11" "E11" "F11" "G11" "H11"
[89] "A12" "B12" "C12" "D12" "E12" "F12" "G12" "H12"
```

```r
summary(data)
```

```
      A01             B01             C01             D01      
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :   0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:   0  
 Median :    0   Median :    0   Median :    0   Median :   0  
 Mean   :   13   Mean   :    6   Mean   :   14   Mean   :   5  
 3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:   2  
 Max.   :34092   Max.   :14704   Max.   :14281   Max.   :3403  
      E01             F01            G01             H01       
 Min.   :    0   Min.   :   0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:   0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :   0   Median :    0   Median :    0  
 Mean   :   14   Mean   :   6   Mean   :    7   Mean   :    4  
 3rd Qu.:    0   3rd Qu.:   0   3rd Qu.:    0   3rd Qu.:    0  
 Max.   :22224   Max.   :5693   Max.   :20038   Max.   :14430  
      A02             B02             C02             D02       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    0  
 Mean   :   14   Mean   :   36   Mean   :   24   Mean   :   16  
 3rd Qu.:    0   3rd Qu.:    4   3rd Qu.:    0   3rd Qu.:    0  
 Max.   :23150   Max.   :68526   Max.   :20262   Max.   :38653  
      E02             F02             G02             H02       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    0  
 Mean   :   12   Mean   :   18   Mean   :   15   Mean   :   15  
 3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    0  
 Max.   :20906   Max.   :38527   Max.   :24591   Max.   :31903  
      A03              B03             C03             D03       
 Min.   :   0.0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:   0.0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :   0.0   Median :    0   Median :    0   Median :    0  
 Mean   :   4.1   Mean   :    9   Mean   :   17   Mean   :   25  
 3rd Qu.:   0.0   3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    1  
 Max.   :2926.0   Max.   :18574   Max.   :33019   Max.   :22298  
      E03             F03             G03             H03       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    0  
 Mean   :   47   Mean   :   30   Mean   :   32   Mean   :   20  
 3rd Qu.:    3   3rd Qu.:    8   3rd Qu.:    0   3rd Qu.:    0  
 Max.   :81906   Max.   :29550   Max.   :60596   Max.   :51778  
      A04              B04             C04             D04      
 Min.   :     0   Min.   :    0   Min.   :    0   Min.   :   0  
 1st Qu.:     0   1st Qu.:    0   1st Qu.:    0   1st Qu.:   0  
 Median :     0   Median :    0   Median :    0   Median :   0  
 Mean   :    43   Mean   :   10   Mean   :   35   Mean   :   3  
 3rd Qu.:     0   3rd Qu.:    1   3rd Qu.:    1   3rd Qu.:   0  
 Max.   :109122   Max.   :11709   Max.   :58601   Max.   :6880  
      E04             F04             G04             H04      
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :   0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:   0  
 Median :    0   Median :    0   Median :    0   Median :   0  
 Mean   :   16   Mean   :    5   Mean   :   14   Mean   :   2  
 3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:   0  
 Max.   :39610   Max.   :19267   Max.   :29707   Max.   :9797  
      A05              B05             C05             D05       
 Min.   :     0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:     0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :     0   Median :    0   Median :    0   Median :    0  
 Mean   :    61   Mean   :   18   Mean   :   20   Mean   :   19  
 3rd Qu.:     1   3rd Qu.:    1   3rd Qu.:    0   3rd Qu.:    0  
 Max.   :189693   Max.   :35899   Max.   :27738   Max.   :22140  
      E05             F05             G05             H05       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    0  
 Mean   :   12   Mean   :   26   Mean   :   28   Mean   :    7  
 3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    0  
 Max.   :45665   Max.   :53400   Max.   :78558   Max.   :19841  
      A06             B06             C06             D06       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    0  
 Mean   :   35   Mean   :    9   Mean   :    8   Mean   :    5  
 3rd Qu.:    2   3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:    0  
 Max.   :72319   Max.   :14724   Max.   :14560   Max.   :11773  
      E06             F06            G06             H06       
 Min.   :    0   Min.   :   0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:   0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :   0   Median :    0   Median :    0  
 Mean   :   21   Mean   :   1   Mean   :   53   Mean   :   18  
 3rd Qu.:    0   3rd Qu.:   0   3rd Qu.:    6   3rd Qu.:    0  
 Max.   :38185   Max.   :5845   Max.   :46241   Max.   :49714  
      A07             B07            C07             D07       
 Min.   :    0   Min.   :   0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:   0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :   0   Median :    5   Median :    0  
 Mean   :   75   Mean   :  13   Mean   :   92   Mean   :   36  
 3rd Qu.:   41   3rd Qu.:   1   3rd Qu.:   55   3rd Qu.:    6  
 Max.   :25005   Max.   :5350   Max.   :33347   Max.   :13164  
      E07              F07             G07             H07        
 Min.   :   0.0   Min.   :    0   Min.   :    0   Min.   :     0  
 1st Qu.:   0.0   1st Qu.:    0   1st Qu.:    0   1st Qu.:     0  
 Median :   0.0   Median :    0   Median :    0   Median :     0  
 Mean   :   5.1   Mean   :   63   Mean   :   63   Mean   :    95  
 3rd Qu.:   0.0   3rd Qu.:   27   3rd Qu.:   30   3rd Qu.:    45  
 Max.   :2407.0   Max.   :25107   Max.   :19423   Max.   :220439  
      A08             B08             C08             D08       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :   10   Median :    0   Median :    0  
 Mean   :   64   Mean   :  129   Mean   :   54   Mean   :   77  
 3rd Qu.:   33   3rd Qu.:   84   3rd Qu.:   26   3rd Qu.:   37  
 Max.   :20502   Max.   :78777   Max.   :24045   Max.   :27117  
      E08             F08            G08             H08       
 Min.   :    0   Min.   :   0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:   0   1st Qu.:    0   1st Qu.:    0  
 Median :    2   Median :   0   Median :    0   Median :    0  
 Mean   :  149   Mean   :  13   Mean   :  104   Mean   :   74  
 3rd Qu.:   70   3rd Qu.:   0   3rd Qu.:   53   3rd Qu.:   39  
 Max.   :92058   Max.   :5133   Max.   :30620   Max.   :37224  
      A09             B09            C09             D09       
 Min.   :    0   Min.   :   0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:   0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :   0   Median :    0   Median :    0  
 Mean   :   95   Mean   :  12   Mean   :   77   Mean   :   72  
 3rd Qu.:   41   3rd Qu.:   1   3rd Qu.:   43   3rd Qu.:   24  
 Max.   :39503   Max.   :4775   Max.   :36453   Max.   :60006  
      E09             F09             G09             H09       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    0  
 Mean   :   98   Mean   :   82   Mean   :   50   Mean   :   66  
 3rd Qu.:   51   3rd Qu.:   32   3rd Qu.:   26   3rd Qu.:   32  
 Max.   :34915   Max.   :33790   Max.   :19967   Max.   :29404  
      A10             B10             C10             D10       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    0  
 Mean   :   34   Mean   :   37   Mean   :   50   Mean   :   88  
 3rd Qu.:    4   3rd Qu.:   16   3rd Qu.:   34   3rd Qu.:   46  
 Max.   :24160   Max.   :19398   Max.   :12113   Max.   :21232  
      E10             F10            G10            H10       
 Min.   :    0   Min.   :   0   Min.   :   0   Min.   :    0  
 1st Qu.:    0   1st Qu.:   0   1st Qu.:   0   1st Qu.:    0  
 Median :    0   Median :   0   Median :   0   Median :    0  
 Mean   :   69   Mean   :  35   Mean   :  10   Mean   :   14  
 3rd Qu.:   35   3rd Qu.:  11   3rd Qu.:   0   3rd Qu.:    0  
 Max.   :23563   Max.   :8293   Max.   :3781   Max.   :13335  
      A11             B11            C11             D11       
 Min.   :    0   Min.   :   0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:   0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :   0   Median :    0   Median :    0  
 Mean   :   82   Mean   :  17   Mean   :   21   Mean   :   43  
 3rd Qu.:   49   3rd Qu.:   0   3rd Qu.:    0   3rd Qu.:   12  
 Max.   :31329   Max.   :4514   Max.   :10313   Max.   :19993  
      E11            F11             G11             H11       
 Min.   :   0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:   0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :   0   Median :    0   Median :    1   Median :    0  
 Mean   :  13   Mean   :   31   Mean   :   87   Mean   :   46  
 3rd Qu.:   0   3rd Qu.:    8   3rd Qu.:   49   3rd Qu.:   20  
 Max.   :7193   Max.   :10857   Max.   :40182   Max.   :33261  
      A12             B12             C12             D12       
 Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:    0   1st Qu.:    0  
 Median :    0   Median :    0   Median :    0   Median :    7  
 Mean   :   88   Mean   :   93   Mean   :   63   Mean   :  119  
 3rd Qu.:   49   3rd Qu.:   52   3rd Qu.:   31   3rd Qu.:   56  
 Max.   :92312   Max.   :27730   Max.   :25082   Max.   :91799  
      E12             F12             G12            H12      
 Min.   :    0   Min.   :    0   Min.   :   0   Min.   :   0  
 1st Qu.:    0   1st Qu.:    0   1st Qu.:   0   1st Qu.:   0  
 Median :    0   Median :    0   Median :   0   Median :   0  
 Mean   :    3   Mean   :    3   Mean   :   2   Mean   :   1  
 3rd Qu.:    0   3rd Qu.:    0   3rd Qu.:   0   3rd Qu.:   0  
 Max.   :10879   Max.   :10650   Max.   :7203   Max.   :4903  
```

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

```
      Chr            Pos           Strand       TrLen         
 2      :1115   Min.   :0.00e+00    : 932   Min.   :      52  
 11     :1112   1st Qu.:3.60e+07   -:6890   1st Qu.:    1833  
 7      :1068   Median :7.45e+07   +:7091   Median :    2962  
        : 932   Mean   :7.58e+07            Mean   :   78581  
 4      : 880   3rd Qu.:1.10e+08            3rd Qu.:    4682  
 5      : 854   Max.   :1.97e+08            Max.   :71445447  
 (Other):8952   NA's   :932                                   
  MinExonHits       MaxExonHits     
 Min.   :      0   Min.   :      1  
 1st Qu.:    119   1st Qu.:    144  
 Median :    599   Median :    656  
 Mean   :   3362   Mean   :   3267  
 3rd Qu.:   1856   3rd Qu.:   1947  
 Max.   :2020769   Max.   :2020771  
                   NA's   :932      
```

Let's characterize this data.
===

- Is the data normal?
- Sparsity / zero-inflation
- Overdispersed

Normality?
===


```r
feature.sum = apply( data, 1, sum )
feature.sum.sorted = sort( feature.sum )
plot( log( feature.sum.sorted ), main = "Total gene counts throughout samples", ylab = "Log total counts", xlab = "Index after sorting" )
abline( v = c(1, 3729,7457,11180,14913), col = c("red", "violet","cyan","blue","green"))
```

![plot of chunk unnamed-chunk-4](scAnalysis-figure/unnamed-chunk-4.png) 

```r
# Min 1, 1st quartile 3729, Median 7457, 3rd quartile 11180, Max 14913
```

Normality?
===


```r
feature.sum.order <- order( feature.sum )
index_min <- feature.sum.order[ 1:10 ]
index_q1 <- feature.sum.order[ 3724:3723 ]
index_median <- feature.sum.order[ 7452:7461 ]
index_q3 <- feature.sum.order[ 11175:11184 ]
index_max <- feature.sum.order[ 14908:14913 ]
plot( density( as.matrix( data[ index_min[1], ] ) ), col = "red" )
```

![plot of chunk unnamed-chunk-5](scAnalysis-figure/unnamed-chunk-5.png) 

Zooming in to Gene Distributions
===
![plot of chunk unnamed-chunk-6](scAnalysis-figure/unnamed-chunk-6.png) 

Sparsity
===


```r
# Check samples for read depth
sample.depth <- apply( data, 2, sum )
summary(sample.depth )
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  10100  183000  337000  556000  941000 2220000 
```
---

```r
barplot( sort( sample.depth ), main = "Sample Depth ", xlab = "Sample", ylab = "Depth")
```

![plot of chunk unnamed-chunk-8](scAnalysis-figure/unnamed-chunk-8.png) 

Sparsity
===


```r
# Percent zero by feature expression
feature.percent.zero <- apply( data, 1, function(x){ length( which(x==0))/length(x)})
feature.mean.no.zero <- apply( data, 1, function(x){ mean( x[ which(x!=0) ] ) })
plot( feature.percent.zero, log( feature.mean.no.zero ), xlab = "log mean", ylab = "Percent Zero", main = "Feature Sparsity by Expression" )
```

![plot of chunk unnamed-chunk-9](scAnalysis-figure/unnamed-chunk-9.png) 

Overdispersion
===

```r
# SD vs Mean
feature.sd.with.zero <- apply( data, 1, sd )
feature.mean.with.zero <- apply( data, 1, mean )
feature.sd.no.zero <- apply( data, 1, function(x){ sd( x[ which(x !=0 )])})
plot( log( feature.mean.no.zero ), log( feature.sd.no.zero ), xlab = "Mean (Log)", ylab = "SD (Log)", main = "SD vs mean (ignoring zeros)" )
```

![plot of chunk unnamed-chunk-10](scAnalysis-figure/unnamed-chunk-101.png) 

```r
plot( log( feature.mean.with.zero ), log( feature.sd.with.zero ), xlab = "Mean (Log)", ylab = "SD (Log)", main = "SD vs mean (ignoring zeros)" )
```

![plot of chunk unnamed-chunk-10](scAnalysis-figure/unnamed-chunk-102.png) 

Can QC Help?
===

- Removing very sparse features
- Imputing outliers

Removing Sparse Features
===


```r
sample.percentile = apply( data, 2, function(x){ quantile(x[x !=0 ], .5)})
feature.noise = which(apply( data, 1, func_min_occurence_at_min_value, sample.percentile ) <= 10)
feature.noise.by.expression = order( apply( data[ feature.noise, ], 1, sum ), decreasing = TRUE)
# What am I removing
plot(density( as.matrix(data[ feature.noise[ feature.noise.by.expression[ 1 ]], ])))
```

![plot of chunk unnamed-chunk-11](scAnalysis-figure/unnamed-chunk-11.png) 

```r
# TODO make a multiple density plot
```

Removing Sparse Features
===


```r
dim( data )
```

```
[1] 14913    96
```

```r
dim( metadata )
```

```
[1] 14913     6
```

```r
data = data[ -1 * feature.noise, ]
metadata = metadata[ -1 * feature.noise ]
dim( data )
```

```
[1] 5654   96
```

```r
dim( metadata )
```

```
[1] 14913     6
```

Sample Read Depth: revisited
===


```r
sample.depth <- apply( data, 2, sum )
depth.colors = polychromatic_palette( length( sample.depth ) )
sample.depth
```

```
    A01     B01     C01     D01     E01     F01     G01     H01     A02 
 187554   93659  202401   65410  192800   91399   99019   52991  201274 
    B02     C02     D02     E02     F02     G02     H02     A03     B03 
 526507  342569  232204  177597  252770  215439  222060   56420  129296 
    C03     D03     E03     F03     G03     H03     A04     B04     C04 
 245420  352991  666297  417858  464753  283792  614285  139734  510160 
    D04     E04     F04     G04     H04     A05     B05     C05     D05 
  41746  230661   76099  200566   29037  873041  249177  282084  275886 
    E05     F05     G05     H05     A06     B06     C06     D06     E06 
 176234  379991  398765  102341  502027  132028  113541   66583  299967 
    F06     G06     H06     A07     B07     C07     D07     E07     F07 
  19785  742681  252405 1041355  181184 1271565  500164   71281  883816 
    G07     H07     A08     B08     C08     D08     E08     F08     G08 
 877940 1316561  890892 1754609  747515 1076466 2084631  182158 1448336 
    H08     A09     B09     C09     D09     E09     F09     G09     H09 
1030883 1332581  160198 1069988 1012292 1367643 1161801  698682  914147 
    A10     B10     C10     D10     E10     F10     G10     H10     A11 
 484067  515619  680259 1226040  965440  472539  148096  187972 1126487 
    B11     C11     D11     E11     F11     G11     H11     A12     B12 
 240928  283019  608197  182873  429258 1212300  636418 1191180 1286476 
    C12     D12     E12     F12     G12     H12 
 878098 1671428   49991   47838   24858   10043 
```
---

```r
barplot( sort(apply( data, 2, sum )))
```

![plot of chunk unnamed-chunk-14](scAnalysis-figure/unnamed-chunk-14.png) 

Can Transforms / Normalization Help?
===

Challenge
===

- We have a large range of read depth per sample
  - Is this a problem?

Waste Not, Want Not
===

- Rarefy?

Dimensionality Reduction and Ordination
===

- A frequently used method is PCA

PCA: in quick theory
===

- Describe PCA

PCA: in practice
===

- Things to be aware of

PCA: in code
===


```r
# Row center and log
data.scaled = t( scale( t( as.matrix( log( data + 1 ) ) ), center=TRUE, scale=TRUE ) )
# Remove constant rows
#TODO data.scaled = data.scaled[, -1 * nearZeroVar( data.scaled ) ]
# Perfrom PCA
results.pca = prcomp( data.scaled, retx = TRUE )
```

PCA: in code
===


```r
plot( results.pca$rotation[,1], results.pca$rotation[,2], pch=16 )
```

![plot of chunk unnamed-chunk-16](scAnalysis-figure/unnamed-chunk-16.png) 

PCA: by depth
===


```r
plot( results.pca$rotation[,1], results.pca$rotation[,2], pch=16, col= depth.colors, xlab="PC1", ylab="PC2", main="PCA by Sample Depth" )
legend( "topright", c(paste(max.pc1," (Max)"), paste(median.pc1," (Median)"), paste(min.pc1," (Min)")), fill=c(depth.colors[max.index], depth.colors[median.index], depth.colors[min.index]) )
```

![plot of chunk unnamed-chunk-18](scAnalysis-figure/unnamed-chunk-18.png) 

Alternatives?
===

- Principle Coordinates Analysis
- Canonical Correspondance Analysis

PCoA: in quick theory
===

PCoA: in practice
===

![Wizard Corgi](images/wizard_corgi.jpg)
- The magic is in the metric
- By default you only get 2 dimensions

PCoA: in code (Bray-curtis)
===


```
Run 0 stress 0.09202 
Run 1 stress 0.09321 
Run 2 stress 0.09257 
Run 3 stress 0.09531 
Run 4 stress 0.1011 
Run 5 stress 0.09594 
Run 6 stress 0.09417 
Run 7 stress 0.09574 
Run 8 stress 0.1986 
Run 9 stress 0.09295 
Run 10 stress 0.09438 
```

```r
plot( nmds.b.c.result$points[,1], nmds.b.c.result$points[,2], col=depth.colors, main="Ordination by Bray Curtis")
legend( "topleft", c(paste(max.pc1," (Max)"), paste(median.pc1," (Median)"), paste(min.pc1," (Min)")), fill=c(depth.colors[max.index], depth.colors[median.index], depth.colors[min.index]) )
```

![plot of chunk unnamed-chunk-20](scAnalysis-figure/unnamed-chunk-20.png) 

PCoA: in code (Jaccard)
===


```r
nmds.j.result = metaMDS( comm=t(pcoa.data), distance="jaccard", k=2, autotransfer=FALSE, trymax=10)
```

```
Run 0 stress 0.09331 
Run 1 stress 0.0956 
Run 2 stress 0.09517 
Run 3 stress 0.09718 
Run 4 stress 0.09445 
Run 5 stress 0.09493 
Run 6 stress 0.0957 
Run 7 stress 0.09435 
Run 8 stress 0.09397 
Run 9 stress 0.09714 
Run 10 stress 0.09428 
```

```r
plot( nmds.j.result$points[,1], nmds.j.result$points[,2], col=depth.colors, main="Ordination by Jaccard")
legend( "topleft", c(paste(max.pc1," (Max)"), paste(median.pc1," (Median)"), paste(min.pc1," (Min)")), fill=c(depth.colors[max.index], depth.colors[median.index], depth.colors[min.index]) )
```

![plot of chunk unnamed-chunk-21](scAnalysis-figure/unnamed-chunk-21.png) 

PCoA: in code (Reciprocal)
===

TODAY

CCA: in quick theory
===

CCA: in practice
===

CCA: in code
===

Next Steps in Ordination
===

Other options to explore
- Sparse PCA

Unsupervised Substructure Discovery
===

Often a goal of scProjects is to describe new structure to a group of cells:
- Heterogeniety of tumor populations
- Novel steps in development
- Robust / dynamic cellular signalling

PCA + Anova: PCA
===


```r
results.pca.features <- prcomp( scale( log( t(data) + 1 ), center=TRUE, scale=TRUE), retx=TRUE)
plot( results.pca.features)
```

![plot of chunk unnamed-chunk-22](scAnalysis-figure/unnamed-chunk-22.png) 

PCA + Anova: PCA
===


```r
hist( abs( results.pca.features$rotation[,1] ), main="Feature Loadings (Abs)", xlab="Feature loadings (Abs)" )
ordered.pca.loadings <- order( results.pca.features$rotation[,1],decreasing = TRUE )
extreme.ends.pca <- c(ordered.pca.loadings[1:250], ordered.pca.loadings[5154:5654])
abline( v=results.pca.features$rotation[,1][extreme.ends.pca], col="#ff000010")
```

![plot of chunk unnamed-chunk-23](scAnalysis-figure/unnamed-chunk-23.png) 

PCA + Anova: visualize
===


```r
data.scaled.subset = t( scale( t( as.matrix( log( data[ extreme.ends.pca, ] + 1 ) ) ), center=TRUE, scale=TRUE ) )
func_heatmap( data.scaled.subset, str_cor ="euclidean")
```

![plot of chunk unnamed-chunk-24](scAnalysis-figure/unnamed-chunk-24.png) 

PCA + Anova: select sample groups
===


```r
#dist.col = dist( t( data.scaled.subset ), method="euclidean")
#dendrogram.col = as.dendrogram( hclust( dist.col, method="average" ) )
gmd.dist <- gmdm2dist( gmdm( t( data.scaled.subset ) ) )
gmd.clust <- css.hclust( gmd.dist, hclust( gmd.dist ) )
gmd.elbow.groups <- elbow.batch( gmd.clust, ev.thres=.9, inc.thres=.05 )
```

PCA + Anova: another view
===


```r
results.pca.subset = prcomp( data.scaled.subset, retx = TRUE )
plot( results.pca.subset$rotation[,1], results.pca.subset$rotation[,2], pch=16, xlab="PC1", ylab="PC2", main="PCA by Sample Depth", col = depth.colors )
```

![plot of chunk unnamed-chunk-26](scAnalysis-figure/unnamed-chunk-26.png) 

PCA + Anova: select feature groups
===

mclust
===

![plot of chunk unnamed-chunk-27](scAnalysis-figure/unnamed-chunk-271.png) ![plot of chunk unnamed-chunk-27](scAnalysis-figure/unnamed-chunk-272.png) ![plot of chunk unnamed-chunk-27](scAnalysis-figure/unnamed-chunk-273.png) ![plot of chunk unnamed-chunk-27](scAnalysis-figure/unnamed-chunk-274.png) ![plot of chunk unnamed-chunk-27](scAnalysis-figure/unnamed-chunk-275.png) 

Compare Methods
===



scOpportunities: cell cycle plot
===

TODAY

scOppotunities: filter by apoptosis
===

TODAY

Statistical Inference in scData
===

- Baysian analysis
- Mixture Models
- Bimodal assumptions

SCDE: in quick theory
===

SCDE: in practice
===

SCDE: in code
===

TODAY

Summary: of the data
===

We are still understanding scData and how to apply it

- Not normal
- Zero-inflated
- Multimodal
- Over-dispersed
...

Summary: of methods
===


Summary: of today
===

- Created expectations on scData
- Explored filtering / normalization techniques
- Applied 3 ordination techniques
- Tried 2 methods to detect substructure
- Applied 1 statistical inference method

Thank you
===

- Aviv Regev
- Alex
- Rahul
- Manik Kuchroo

Questions?
===

![Answers Corgi](images/graduate_corgi.jpg)

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
