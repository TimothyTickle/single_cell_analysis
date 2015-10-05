###########################################
# Data QC and Normalization methods
###########################################

filter.by.occurence <- function(
### Filter genes (rows) by require those kept to have atleast a min values in a min number of samples.
### This avoids using averages and it stable when adding samples to the study.
df_data,
### Data frame to be filtered
d_min_value,
### Minimum value a measurement of gene must be to be counted in a sample.
d_min_occurence
### Minimum times a gene must be found in a sample (at the d_min_value) in order to not be filtered.
){
  vctr_f_keep = apply( df_data, 1, function( x ){ return( length( which( x >= d_min_value ) ) >= d_min_occurence ) } )
  return( df_data[ vctr_f_keep,] )
}


normalize.cpx <- function(
### Normalize counts to counts per a X where X is a median depth of the sample.
### Normalizes within columns
df_data,
### Count data to be transformed
...
){
  d_median_sums = median( apply( df_data, 2, sum ) )
  d_magnitude = as.integer(paste( c(c( "1"), rep( "0", nchar(as.character(d_median_sums)))),collapse=""))
  return( log2( ( normalize.tss( df_data ) * d_magnitude ) + 1 ))
}


normalize.tss <- function(
### Total Sum Scaled
### Normalize columns of the data set buy dividing each observation by the total of the sample
df_data,
### Data to be transformed
...
){
  return( sweep( df_data, 2, colSums( df_data), "/" ) )
}


plot.saturation.curve <- function(
### in Silico bootstrap cell complexity at varying levels to get a saturation curve for sample depth
vctr_values,
### Sample expression
i_depth_increment,
### Depth increment for sampling from 0 to d_max_depth 
d_max_depth = NA,
### Make depth to explore (if NA will be set to the depth of the sample).
i_iterations = 10
### How many times a specific depth is sampled (increasing this increases the stability of each depth's measurement on the graph). 
){
  # Set up a population to sample from
  vctr_i_population = c()
  for( i_values in 1:length(vctr_values) ){
    vctr_i_population = c(vctr_i_population, rep(i_values, vctr_values[i_values] ) )
  }
  vctr_bootstrapped_genes_expressed = c()
  
  # If max depth is not given, use the sum of the values
  if( is.na( d_max_depth )){
    d_max_depth = as.integer(sum( vctr_values ))
  }
  
  # Increment can not be bigger than the max depth given
  i_depth_increment = min(i_depth_increment, d_max_depth)
  
  # Set sequencing of depths sampled
  vctr_d_depths = seq( i_depth_increment, d_max_depth, i_depth_increment )
  
  # For each depth increment
  for( d_depth in vctr_d_depths )
  {
    vctr_exprs_genes = c()
    for( i_iter in 1:i_iterations )
    {
      # Bootstrap the counts
      # Count how many genes express
      vctr_exprs_genes = c( length( unique( sample( vctr_i_population, d_depth, replace = TRUE ))), vctr_exprs_genes )
    }
    vctr_bootstrapped_genes_expressed = c( vctr_bootstrapped_genes_expressed, mean(vctr_exprs_genes) )
  }
  # Add Zero, zero to anchor plot
  vctr_d_depths=c(0,vctr_d_depths)
  vctr_bootstrapped_genes_expressed=c(0,vctr_bootstrapped_genes_expressed)
  # Plot
  plot( vctr_d_depths, vctr_bootstrapped_genes_expressed, 
        main="Saturation Curve", xlab="Depth (in Silico)", ylab="Genes Detected",
        pch=16, col="#0000ff55")
}


plot.cell.complexity <- function(
  ### This is only a demonstration of how outlier samples can be identified by complexity.
  ### Plot cell complexity and identify outliers
  vctr_values,
  ### Cell complexity
  ...
){
  bx.qt = as.integer( boxplot( vctr_values, plot=FALSE )$stats )
  lower.outlier = which( vctr_values < bx.qt[1] )
  higher.outlier = which( vctr_values > bx.qt[5] )
  vioplot( vctr_values, rectCol="white", col="cyan", colMed="black", border="purple" )
  stripchart( vctr_values, add=TRUE, vertical=TRUE, pch=21, jitter= .1, method="jitter")
  title( main="Cell Complexity", xlab="Study", ylab="Cell Complexity" )
  return( c( lower.outlier, higher.outlier ) )
}


demo.complexity.outliers <- function(
  ### Returns a list of values that demostrate low outliers in complexity.
  ### This is just for the presentation.
){
  return(c(200,300,212,354,213,2252,1633,1667,2162,2674,1819,1159,3176,2902,2827,2349,3001,1735,1904,1354,1952,1775,1799,1704,2410,2578,1935,2770,1914,1776,1428,2380,1773,2965,1978,1539,1371,1305,563,1583,1135,3164,1661,2300,1609,1997,2361,2750,1386,1787,1536,1577,2236,2681,2071,2278,1768,1578,2342,1775,2492,3070,4451,1785,1943,2084,1833,2969,1918,2882,2236,2946,1673,1407,4029,2731,2615,2262,2674,2081,1179,2301,2354,1885,1524,2766,1858,2153,1849,1387,1710,2137,2047,1738,2814,2131,1075,2187,2058,2667,2410,2125,2228,3330,2696))
}


#####################################
# Seurat Associated
#####################################

read.into.seurat <- function(
### Read data into a seurat object
str_data_file,
### Text file of expression
...
){
  # Read table
  nbt.data=read.table( str_data_file, ... )
  # Log and add 1
  nbt.data=log( nbt.data + 1 )
  # Make Seurat data file
  return( new( "seurat", raw.data=nbt.data ) )
}


prep.pca.seurat <- function(
  ### Wrapper function to simpily prep for PCA.
  data.seurat,
  ### Data to prep
  y.cutoff,
  ### Dispersion lowest value cut-off
  x.low.cutoff
  ### Average expression lowest value cut-off
){
  data.seurat=mean.var.plot( data.seurat, y.cutoff=y.cutoff, x.low.cutoff=x.low.cutoff, fxn.x=expMean, fxn.y=logVarDivMean, do.plot=FALSE)
  data.seurat=pca( data.seurat, do.print=FALSE)
  return( data.seurat )
}
