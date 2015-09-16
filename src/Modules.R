
# Constants
C_STR_DENSITY_COLOR = "cyan" #"darkmagenta" # #BB00BB
### Color
STR_COLOR_AFTER = "#00660075" # green
### Used for after a method (used as a secondary to STR_COLOR_BEFORE)
STR_COLOR_BEFORE = "#DD640075" # orange
### USed for the before color in a before / after picture.
### Also used as the selection color for plots.
STR_COLOR_STANDARD = "#00000075" # transparent black ( grey )
### Color to be used for normal plotting


################################
# Color palettes
################################

monochromatic_palette <- function(
### Return the standard monochromatic color scheme
i_number_colors = 100
){
  return( colorpanel( i_number_colors, "white", "violet", "purple" ) )
}

polychromatic_palette <- function(
### Return the standard polychromatic color scheme
i_number_colors = 100
){
  return( colorpanel( i_number_colors, "purple", "black", "yellow" ) )
}

metadata_palette <- function(
### Return the standard color scheme for metadata
i_number_colors = 100
){
  return( rainbow( i_number_colors ) )
}

min_occurence_at_min_value <- function(
### Returns if a vector has a minimum occurence of values at or equal to the given min value
### This is a helper function for func_filter_by_percentile_percentage
vctr_values_to_evaluate,
### Values to check if the are greater or equal to the min value
vctr_min_values,
### Minimum values (positionally) a value must be to be counted
d_min_occurence,
### Minimum counts a feature must have above the value
...
){
  return( length( which( vctr_values_to_evaluate >= vctr_min_values ) ) )
}

############################################
# Data plots
############################################
plot.quantiles <- function(
  ### Plot data arbitrarly selected for the quantiles of the data
  data
){
  # Get the sorted order for the genes
  feature.sum.order <- order( apply( data, 1, sum) )
  feature.summary <- summary( 1:length(feature.sum.order))
  # Get the 10 most sparse genes
  index_min <- feature.sum.order[ 1:10 ]
  # 1 quartile sparsity genes
  index_q1 <- feature.sum.order[ floor(feature.summary[[2]]):(floor(feature.summary[[2]])+9) ]
  # Median sparsity genes
  index_median <- feature.sum.order[ floor(feature.summary[[3]]):(floor(feature.summary[[3]])+9) ]
  # 3rd quartile sparsity genes
  index_q3 <- feature.sum.order[ floor(feature.summary[[5]]):(floor(feature.summary[[5]])+9) ]
  # Get the least sparse genes
  index_max <- feature.sum.order[ (length(feature.sum.order)-9):length(feature.sum.order) ]
  plot( x=0,y=0,type="p", xlim=c(0,log(max(data)+1)), ylim=c(0,.75), main="Gene Count Distributions by Sparsity", ylab="Density of gene counts", xlab="Gene count value (Log)" )
  for( i_q3_plot in index_q3 ){ lines( density( as.matrix( log(data[ i_q3_plot, ]+1))), col = "#0000ff75")}
  for( i_median_plot in index_median ){ lines( density( as.matrix( log(data[ i_median_plot, ]+1))), col = "#00ffff75")}
  for( i_min_plot in index_min ){ lines( density( as.matrix( log(data[ i_min_plot, ]+1))), col = "#ff000075")}
  for( i_q1_plot in index_q1 ){ lines( density( as.matrix( log(data[ i_q1_plot, ]+1))), col = "#EE82EE75")}
  for( i_max_plot in index_max ){ lines( density( as.matrix( log(data[ i_max_plot, ]+1))), col = "#00ff0075")}
  legend( "topright", c("Min","Q1","Median","Q3","MAX"), fill=c("#ff0000","#EE82EE","#00ffff","#0000ff","#00ff00"), title="Sparsity group" )
  
}

###########################################
# Data QC and Normalization methods
###########################################

func_filter_by_occurence <- function(
df_data,
d_min_value,
d_min_occurence
){
  vctr_f_keep = apply( df_data, 1, function( x ){ return( length( which( x >= d_min_value ) ) >= d_min_occurence ) } )
  return( df_data[ vctr_f_keep,] )
}

func_cpx <- function(
### Normalize counts to counts per million
### Normalizes within columns
df_data,
### Count data to be transformed
...
){
  d_median_medians = median( apply( df_data, 2, median ) )
  d_magnitude = as.integer(paste( c(c( "1"), rep( "0", nchar(as.character(d_median_medians)))),collapse=""))
  return( log2( ( func_tss( df_data ) * d_magnitude ) + 1 ))
}


func_tss <- function(
### Total Sum Scaled
### Normalize columns of the data set buy dividing each observation by the total of the sample
df_data,
### Data to be transformed
...
){
  return( sweep( df_data, 2, colSums( df_data), "/" ) )
}

func_plot_saturation_curve <- function( 
vctr_values,
i_depth_increment,
d_max_depth = NA,
i_iterations = 10
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
        main="Saturation Curve", xlab="Depth", ylab="Genes Detected",
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
  lower.outlier = which( vctr_values <= bx.qt[1] )
  higher.outlier = which( vctr_values >= bx.qt[5] )
  vioplot( genes.per.cell, rectCol="white", col="cyan", colMed="black", border="purple" )
  stripchart( genes.per.cell, add=TRUE, vertical=TRUE, pch=21, jitter= .1, method="jitter")
  title( main="Cell Complexity", xlab="Study", ylab="Cell Complexity" )
  return( c( lower.outlier, higher.outlier ) )
}

###########################################
# Unsupervised ordination
###########################################


heatmap <- function(
### Create a heatmap of data frame
df_frame,
### Data for which to show a heatmap
str_cor = "euclidean",
### The correlation metric to use
str_linkage = "average",
### Any valid linkage method for hclust
vctr_grouping = NULL,
### Optional Groups the samples
str_output_dir = NULL,
### The directory to output figure and such
plt_colors = NULL,
### Color pallette, if not given polychromatic default will be used
str_title = "Data Heatmap",
### Title of plot
...
){
  ## TODO Tranform (center and such) should happen before this

  # Get a consistent color scheme for the plotting
  if( is.null( plt_colors ) )
  {
    plt_colors = polychromatic_palette( )
  }

  # Get distance matrix
  dist_row = dist( df_frame, method = str_cor )
  dist_col = dist_row
  if( !nrow( df_frame ) == ncol( df_frame ) )
  {
    dist_col = dist( t( df_frame ), method = str_cor )
  }

  # Visualize heatmap and possibly record in pdf file
  if( ! is.null( str_output_dir ) )
  {
    pdf( paste( str_output_dir, "data_heatmap.pdf" ) )
  }

  if( !is.null( vctr_grouping ) )
  {
    # Set the groupings to colors
    lret_coloring = func_factor_to_metadata_color( vctr_grouping )
    vctr_grouping_colors = matrix( lret_coloring$vctr_grouping_colors, nrow = 1 )
    names( vctr_grouping_colors ) = names( lret_coloring$vctr_grouping_colors )
    heatmap.3b( x=df_frame, dendrogram='both', 
              Rowv=as.dendrogram( hclust( dist_row ), method = str_linkage ), 
              Colv=as.dendrogram( hclust( dist_col ), method = str_linkage ), 
              scale='none', symm=TRUE, key=TRUE, col = plt_colors,
              density.info='histogram', denscol = C_STR_DENSITY_COLOR, trace='none', 
              symkey=FALSE, margins=c(10,10), cexCol=1, ColSideColors = vctr_grouping_colors,
              cexRow=1, cex.main=0.75, main = str_title, ...)
  } else {
    heatmap.3b( x=df_frame, dendrogram='both',
              Rowv=as.dendrogram( hclust( dist_row ), method = str_linkage ),
              Colv=as.dendrogram( hclust( dist_col ), method = str_linkage ),
              scale='none', symm=TRUE, key=TRUE, col = plt_colors,
              density.info='histogram', denscol = C_STR_DENSITY_COLOR, trace='none',
              symkey=FALSE, margins=c(10,10), cexCol=1,
              cexRow=1, cex.main=0.75, main = str_title, ...)
  }

  if( ! is.null( str_output_dir ) )
  {
    dev.off()
  }
}


cor_matrix <- function(
### Show the correlation of the date in the data frame (samples or features)
df_frame,
### Data for which to show correlation
str_cor = "spearman",
### The correlation metric to use
vctr_grouping = NULL,
### Grouping of samples or features by color
str_linkage = "average",
### Any valid linkage method for hclust
str_output_dir = NULL,
### The directory to output figure and such
f_correlate_samples = TRUE,
### Indicates if samples (TRUE; columns) or features (FALSE; rows) are correlated
str_title = NULL,
### Title for plot
...
){
  # Get a consistent color scheme for the plotting
  plt_colors = polychromatic_palette( )

  # Create correlation matrix
  mtrx_cor = NA
  if( f_correlate_samples )
  {
    mtrx_cor = cor( df_frame, method = str_cor, use = 'pairwise.complete.obs' )
  } else {
    mtrx_cor = cor( t( df_frame ), method = str_cor, use = 'pairwise.complete.obs' )
  }

  # Remove diag
  diag( mtrx_cor ) = NA

  # Distance between samples
  dgrm_samples = as.dendrogram( hclust( as.dist( 1-mtrx_cor ), method = str_linkage ) )

  # Visualize heatmap and possibly record in pdf file
  if( ! is.null( str_output_dir ) )
  {
    str_file_name = "sample_correlation_matrix.pdf"
    if( ! f_correlate_samples )
    {
      str_file_name = "feature_corrleation_matrix.pdf"
    }
    pdf( paste( str_output_dir, str_file_name ) )
  }

  if( is.null( str_title ) )
  {
    str_title = "Sample Correlation Matrix"
    if( ! f_correlate_samples )
    {
      str_title = "Feature Correlation Matrix"
    }
  }

  if( !is.null( vctr_grouping ) )
  {
      lret_grouping_colors = func_factor_to_metadata_color( vctr_grouping )
      vctr_grouping_colors = matrix( lret_grouping_colors$vctr_grouping_colors, nrow = 1 )
      names( vctr_grouping_colors ) = names( lret_grouping_colors$vctr_grouping_colors )
      heatmap.3b( mtrx_cor, dendrogram = "both", Rowv = dgrm_samples, Colv = dgrm_samples, col = plt_colors, scale = "none", revC = TRUE, key = TRUE, density.info = "histogram", trace = "none", margin.for.labRow = 10, margin.for.labCol = 10, cexCol = 1, cexRow = 1, cex.main = 0.75, main = str_title, sep.color = "black", sep.lwd = 1, denscol = C_STR_DENSITY_COLOR, ColSideColors = vctr_grouping_colors, RowSideColors = t( vctr_grouping_colors ), ... ) 
  } else {
      heatmap.3b( mtrx_cor, dendrogram = "both", Rowv = dgrm_samples, Colv = dgrm_samples, col = plt_colors, scale = "none", revC = TRUE, key = TRUE, density.info = "histogram", trace = "none", margin.for.labRow = 10, margin.for.labCol = 10, cexCol = 1, cexRow = 1, cex.main = 0.75, main = str_title, sep.color = "black", sep.lwd = 1, denscol = C_STR_DENSITY_COLOR, ... ) 
  }

  if( ! is.null( str_output_dir ) )
  {
    dev.off()
  }
}


func_factor_to_metadata_color <- function(
vctr_factors
### Factors to turn to metadata colors
){
  vctr_grouping_colors = rep("NA", length( vctr_factors ) )
  names( vctr_grouping_colors ) = as.character( vctr_factors )
  vctr_levels = levels( vctr_factors )
  vctr_colors = metadata_palette( length( vctr_levels ) )
  for( i_color_index in 1:length( vctr_levels ) )
  {
    vctr_grouping_colors[ which( vctr_factors == vctr_levels[ i_color_index ] ) ] = vctr_colors[ i_color_index ]
  }
  return( list( vctr_grouping_colors = vctr_grouping_colors, vctr_levels = vctr_levels, vctr_colors = vctr_colors ) )
}

#####################################
# Seurat Associated
#####################################

read.into.seurat <- function(
  str_data_file,
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
  y.cutoff,
  x.low.cutoff
){
  nbt=mean.var.plot(nbt,y.cutoff=y.cutoff,x.low.cutoff=x.low.cutoff,fxn.x = expMean,fxn.y = logVarDivMean,do.plot=FALSE)
  nbt=pca(nbt,do.print=FALSE)
  return( nbt )
}

#####################################
### Monocle associated
#####################################

make_cell_data_set <- function(
expression_file,
cell_phenotype_file,
gene_metadata_file
){
  # Read in expression as matrix
  monocle.exprs <- as.matrix( read.table(expression_file ) )
  # Read in cell phenotype
  monocle.cell.meta <- read.table( cell_phenotype_file )
  # Read in gene metadata
  monocle.gene.meta <- read.table( gene_metadata_file )
  # Prep cell phenotype data
  monocle.pheno <- new( "AnnotationDataFrame", data=monocle.cell.meta )
  # Prep gene metadata
  monocle.feature.data <- new( "AnnotatedDataFrame", data=monocle.gene.meta )
  # Hold everything in a CellDataSet object
  return( newCellDataSet( monocle.exprs, phenoData=monocle.pheno, featureData=monocle.feature.data ) )
}

get_monocle_presentation_data <- function(){
  data( HSMM_expr_matrix )
  data( HSMM_gene_annotation )
  data( HSMM_sample_sheet )
  monocle.pheno <- new( "AnnotatedDataFrame", data=HSMM_sample_sheet )
  monocle.feature.data <- new( "AnnotatedDataFrame", data=HSMM_gene_annotation )
  return( newCellDataSet( as.matrix(HSMM_expr_matrix), phenoData=monocle.pheno, featureData=monocle.feature.data ) )
}

get_monocle_presentation_marker_genes <- function(){
  return( c("MEF2C","MEF2D","MYF5","ANPEP","PDGFRA","MYOG","TPM1","TPM2","MYH2","MYH3","NCAM1","TNNT1","TNNT2","TNNC1","CDK1","CDK2","CCNB1","CNNB2","CCND1","CCNA1","ID1") )
}

plot_log_normal_monocle <- function(
  data.set
){
  # Log transform data
  data.logged <- log(exprs(data.set[ monocle.expr.genes,]))
  # Standardize each gene
  data.logged <- t(scale(t(data.logged)))
  # Plot
  qplot(value, geom="density", data=melt(data.logged)) + stat_function(fun=dnorm, size=0.5, color="red") + xlab("Standarized log") + ylab("Density")
}

select_ordering_genes <- function(
data.set,
genes.expressed,
genes.of.interest,
model.formula,
qvalue.threshold
){
  data.set.subset <- subset_to_genes( data.set, genes.of.interest )
  diff.marker.genes <- differentialGeneTest(data.set.subset, fullModelFormulaStr=model.formula)
  sig.marker.genes <- row.names(subset(diff.marker.genes, qval < qvalue.threshold))
  return( intersect( sig.marker.genes, genes.expressed ) )
}

order_cells_wrapper <- function( 
data.set, 
genes, 
use_irlba, 
num_paths, 
reverse 
){
  data.set <- setOrderingFilter( data.set, genes )
  data.set <- reduceDimension( data.set, use_irlba=use_irlba )
  return( orderCells( data.set, num_paths=num_paths, reverse=reverse) );
}

subset_to_genes <- function(
data.set,
genes
){
  return( data.set[ row.names( subset(fData( data.set ), gene_short_name %in% genes ) ),] )
}
