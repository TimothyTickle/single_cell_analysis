
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


###########################################
# Normalization methods
###########################################


func_cpm <- function(
### Normalize counts to counts per million
### Normalizes within columns
df_data,
### Count data to be transformed
...
){
  return( t( t( df_data )/colSums( df_data ) ) * 1e6 )
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


func_upper_quartile <- function(
### Normalize values using the upper quartile value
### http://vinaykmittal.blogspot.com/2013/10/fpkmrpkm-normalization-caveat-and-upper.html
df_data,
### Count data to be transformed
f_scale = FALSE,
...
){
  loginfo( "Upper quartile normalization" )

  df_normalized = df_data

  # Remove transcripts / genes/ rows that have 0 expression in all samples
  vctr_zero_rows = which( apply( df_data, 1, sum ) == 0 )
  if( length( vctr_zero_rows ) )
  {
    df_normalized = df_normalized[ -1 * vctr_zero_rows, ]
  }

  # For each column / sample find the 75th percentile value (upper quartile)
  vctr_upper_quartile = apply( df_normalized, 2, quantile, probs = c(0.75) ) 

  # Divide all the expression values by the upper quartile
  df_normalized = sweep( df_normalized, MARGIN = 2, vctr_upper_quartile, "/" ) 

  # Optionally scale all the values by the mean of the 75th percentile quartiles
  # To increase very low values created by the normalization
  if( f_scale )
  {
    return( df_normalized * mean( vctr_upper_quartile ) )
  }

  return( df_normalized )
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


#############################
# Distance metrics
#############################

reciprocal_weighting <- function(
x,
### Expected matrix (rows as feature)
groups = NULL,
### Groupings of data
k = 0.95,
### Max value
n.cores = 1
### Number of cores to run on
){
  # Fit error model
  o.ifm <- scde.error.models( counts=x, groups=groups, n.cores=n.cores, threshold.segmentation=TRUE, save.crossfit.plots=FALSE, save.model.plots=FALSE, verbose=1 )
  # Filter out poor fits
  o.ifm <- o.ifm[ o.ifm$corr.a > 0, ]

  # Get expression maginitude estimates
  o.fpm <- scde.expression.magnitude(o.ifm,counts=x)

  # Get cell names to iterate through
  cell.names = names( x )
  require(boot)

  # Make distance matrix
  reciprocal.dist <- as.dist( 1 - do.call( rbind, mclapply( cell.names, function(nam1)
  {
    unlist(lapply(cell.names,function(nam2)
    {
      # Reciprocal probabilities
      f1 <- scde.failure.probability( models=o.ifm[nam1,,drop=FALSE], magnitudes=o.fpm[ ,nam2 ])
      f2 <- scde.failure.probability( models=o.ifm[nam2,,drop=FALSE], magnitudes=o.fpm[ ,nam1 ])
      # Weight factor
      pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
      boot::corr(log10(cbind(data[,nam1],data[,nam2])+1),w=pnf)
    }))
  },mc.cores=n.cores)),upper=FALSE)
  return( reciprocal.dist )
}


#####################################
### Misc
#####################################
stick_breaking <- function(
i_variance_count
### Number of breaks
){
  vctr_break_distribution = c()
  for( i_k_break in 1:i_variance_count )
  {
    vctr_break_distribution = c(vctr_break_distribution, (1/i_variance_count) * sum( 1 / ( i_k_break : i_variance_count ) ) )
  }
  return( vctr_break_distribution )
}
