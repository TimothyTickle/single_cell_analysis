
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


# Key words for modules listed below
vctr_batch_methods = c( "none", "combat", "edgeR", "regression", "PCA_skimming" )
vctr_discriminate_genes_methods = c( "none", "anova" )
vctr_filters_methods = c( "none", "percentile", "sd", "downsample", "head", "occurence", "sparsity" )
vctr_impute_methods = c( "none", "complete_case", "mean", "median", "multiple_imputation" )
vctr_normalization_methods = c( "none", "cpm", "upper_quartile" )
vctr_ordination_methods = c( "none", "pca", "cca", "nmds", "scde" )
vctr_select_sample_groups_methods = c( "none", "anova", "pca" )
vctr_select_feature_group_methods = c( "none", "anova", "pca" )
vctr_select_transform = c( "none", "asin_sqrt", "center", "log", "tss", "zscale" )

# Custom metrics
vctr_custom_distance_metric = c( "reciprocal" )

##################################
# Description of data
##################################


func_describe_metadata <- function(
### Describes the matrix of metadata
df_metadata,
### Output project directory
str_output_dir
){
  loginfo( "Describe metadata" )

  # Manage output directory for PDF
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }
  pdf( file.path( str_output_dir, "describe_metadata.pdf"), useDingbats=FALSE )

  # Counts of metadata
  vctr_missing = apply( df_metadata, 2, function(x) length( union(which(is.na(x)), which(is.nan(x)))))
  vstr_labels = paste( names(df_metadata), "(", nrow( df_metadata ), ")", sep = " ")
  barplot( vctr_missing, names.arg = vstr_labels, ylab = "Count of missing values", xlab = "Metadata", main = "Missing Values per Metadata" )
  dev.off()

  # Metadata associations
  if( sum( vctr_missing ) > 0 )
  {
    i_ncol = ncol( df_metadata )
    i_nrow = nrow( df_metadata )
    df_metadata_discrete = as.data.frame( matrix(rep(0, i_ncol * i_nrow), nrow = i_nrow, ncol = i_ncol ) )
    names( df_metadata_discrete ) = names( df_metadata )
    row.names( df_metadata_discrete ) = row.names( df_metadata )
    df_metadata_discrete[ is.na( df_metadata ) ] = 1
    dist_cor = 1 - vegdist( t(df_metadata_discrete), method = "jaccard" )
    df_metadata_discrete =  data.frame( t(combn( names( df_metadata_discrete ), 2 ) ), as.numeric( dist_cor ))
    func_write_text( "Similarity of Missing Data\n\n", file.path( str_output_dir, "describe_metadata.txt" ))
    func_write_frame( df_metadata_discrete, file.path( str_output_dir, "describe_metadata.txt" ), TRUE)
  }
}


func_describe_count_data_matrix <- function(
### Make basic plots to describe the data
df_frame,
### Data to describe
df_groupings,
### Metadata to describe the data
str_output_dir,
### Directory to output data and visualizations
str_cor_metric = "euclidean",
### Correlation metric to correlate samples
vctr_color_groups,
### Coloring for sample metadata
...
){
  loginfo( "Describe Count Data. Please note assumptions are made here that the data is count." )

  # Manage output directory for PDF
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }
  pdf( file.path( str_output_dir, "describe_count_data.pdf"), useDingbats=FALSE )

  # Basic counts
  plot.new()
  func_describe_data_frame_values( df_frame )

#  # Bulk values
#  vctr_all_values_logged = as.vector( as.matrix( func_log( df_frame ) ) )
#  vctr_all_values = as.vector( as.matrix( df_frame ) )
#  i_all_values_length = length( vctr_all_values )
#  vctr_non_zero = vctr_all_values[ which( vctr_all_values != 0 ) ]
#  i_non_zero_values_length = length( vctr_non_zero )
#
#  # General distribution
#  hist( vctr_all_values, breaks=min( 200, i_all_values_length ), main="General Distribution" )
#  hist( vctr_all_values_logged, breaks=min( 200, i_all_values_length ), main="General Distribution (Logged)" )
#  vctr_all_values_logged = NULL
#
#  # Get distribution without zero values
#  hist( vctr_non_zero, breaks=min( 200, i_non_zero_values_length ), main="General Distribution ( Zeros removed )" )
#  hist( log( vctr_non_zero + 1 ), breaks=min( 200, i_non_zero_values_length ), main="General Distribution (Logged with zeros removed )" )

  # Heat map raw data distribution
  # Add metadata as columns
  func_heatmap( df_frame, str_cor = str_cor_metric, str_linkage = "average", vctr_factor_groups = NULL, plt_colors = func_monochromatic_palette( ),... ) 

  # Heat map correlation matrix
  # TODO diag = NA for copy matrix
  func_cor_matrix( df_frame, str_cor = "spearman", str_linkage = "average", f_correlate_samples = TRUE, ... )
  func_cor_matrix( df_frame, str_cor = "spearman", str_linkage = "average", f_correlate_samples = FALSE, ... )

  # Sample distributions
  boxplot( df_frame, main="Sample distributions" )
  df_frame_logged = func_log( df_frame )
  boxplot( df_frame_logged, main="Logged Sample Distributions" )

  # Plot the means vs the sd
  # Plot the number of zeros by the mean
  vctr_means = apply( df_frame, 2, mean )
  vctr_means_no_zeros = apply( df_frame, 2, function( x ) mean( x[ which( x != 0 ) ] ) )
  vctr_sd = apply( df_frame, 2, sd )
  vctr_sd_no_zeros = apply( df_frame, 2, function(x) sd( x[ which( x != 0 ) ] ) )
  vctr_zeros = apply( df_frame, 2, function( x ) length( which( x == 0 ) ) )
  plot( vctr_means, vctr_sd, main = "Mean vs SD", xlab = "Mean", ylab = "Standard deviation" )
  plot( vctr_means_no_zeros, vctr_sd_no_zeros, main = "Mean vs SD (No Zeros)", xlab = "Mean (No zeros)", ylab = "Standard deviation (No zeros)" )

  # Plot number of zeros vs means excluding zero
  plot( vctr_means_no_zeros, vctr_zeros, main = "Mean (No zeros) vs Number Zeros", xlab = "Mean ( No zeros)", ylab = "Number zeros" )

  # Plot by metadata
  # TODO Use metadata properly
  if( ! is.null( df_groupings ) )
  {
    par( mar = c( 5.1, 4.1, 4.1, 8.1 ), xpd = TRUE )
    for( str_metadata in names( df_groupings ) )
    {
      # Make factor data
      vctr_groups = df_groupings[[ str_metadata ]]
      vctr_grouping_colors = rep( NA, length( vctr_groups ) )
      if( is.numeric( vctr_groups ) || is.integer( vctr_groups ) )
      {
        vctr_groups = cut( vctr_groups, min( 10, length( vctr_groups ) ) )
      } else {
        vctr_groups = factor( vctr_groups )
      }

      # Change to colors
      lret_grouping_colors = func_factor_to_metadata_color( vctr_groups )

      # Add color to the plots
      boxplot( df_frame, main="Sample distributions", col = lret_grouping_colors$vctr_grouping_colors, las = 2 )
      legend( "topright", inset = c( -0.2, 0 ), title = str_metadata, legend = lret_grouping_colors$vctr_levels, col = lret_grouping_colors$vctr_colors, pch = 16 ) 
      boxplot( df_frame_logged, main="Logged Sample Distribution", col = lret_grouping_colors$vctr_grouping_colors, las = 2 )
      legend( "topright", inset = c( -0.2, 0 ), title = str_metadata, legend = lret_grouping_colors$vctr_levels, col = lret_grouping_colors$vctr_colors, pch = 16 )
      plot( vctr_means, vctr_sd, main = "Mean vs SD", xlab = "Mean", ylab = "Standard deviation", pch = 16, col = lret_grouping_colors$vctr_grouping_colors )
      legend( "topright", inset = c( -0.2, 0 ), title = str_metadata, legend = lret_grouping_colors$vctr_levels, col = lret_grouping_colors$vctr_colors, pch = 16 )
      plot( vctr_means_no_zeros, vctr_sd_no_zeros, main = "Mean vs SD (No Zeros)", xlab = "Mean (No zeros)", ylab = "Standard deviation (No zeros)", pch = 16, col = lret_grouping_colors$vctr_grouping_colors )
      legend( "topright", inset = c( -0.2, 0 ), title = str_metadata, legend = lret_grouping_colors$vctr_levels, col = lret_grouping_colors$vctr_colors, pch = 16 )
      plot( vctr_means_no_zeros, vctr_zeros, main = "Mean (No zeros) vs Number Zeros", xlab = "Mean ( No zeros)", ylab = "Number zeros", col = lret_grouping_colors$vctr_colors )
      legend( "topright", inset = c( -0.2, 0 ), title = str_metadata, legend = lret_grouping_colors$vctr_levels, col = lret_grouping_colors$vctr_colors, pch = 16 )

      # Heat map raw data distribution
      # Add metadata as columns
      func_heatmap( df_frame, str_cor = str_cor_metric, str_linkage = "average", plt_colors = func_monochromatic_palette( ), vctr_grouping = vctr_groups, ... ) 
      # Heat map correlation matrix
      func_cor_matrix( df_frame, str_cor = "spearman", str_linkage = "average", f_correlate_samples = TRUE, vctr_grouping = vctr_groups, ... )
      func_cor_matrix( df_frame, str_cor = "spearman", str_linkage = "average", f_correlate_samples = FALSE, ... )
    }
  }
  dev.off()
}


func_describe_data_frame_values = function(
### Prints basic counts and stats about a data frame / matrix
df_frame
### The data frame / matrix to describe
){
  # Calculate basic descriptions
  i_nrows = nrow( df_frame )
  i_cols = ncol( df_frame )
  vctr_row_counts = apply( df_frame, 1, sum )
  vctr_col_counts = apply( df_frame, 2, sum )
  i_zero_rows = length( which( vctr_row_counts == 0 ) )
  i_zero_cols = length( which( vctr_col_counts == 0 ) )

  # Report values
  i_line_increment = 0
  mtext( "Describing count data:", side=3, adj = 0 ,line = i_line_increment )
  i_line_increment = i_line_increment - 2
  str_line = paste( "Number of Rows:", i_nrows )
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
  str_line = paste( "Number of Rows with no measurments:", i_zero_rows )
  i_line_increment = i_line_increment - 1
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
  str_line = paste( "Number of columns:", i_cols )
  i_line_increment = i_line_increment - 1
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
  str_line = paste( "Number of columns with no measurements:", i_zero_cols )
  i_line_increment = i_line_increment - 1
  mtext( str_line, side = 3, adj = 0, line = i_line_increment )
}


######################################
# Data transformation
######################################

func_log <- function(
### Transform the given data frame to log2
### Expects values of ) or greater
### 1 is added to the value before transformation
### Transforms within columns
df_data
### Data to be transformed
){
  return( log2( df_data + 1 ) )
}


func_zscale <- function(
### Mean centers and standardizes by sd to make norm values
### Transforms within columns
df_data
### Data to be transformed
){
  return( data.frame( scale( df_data, center = TRUE, scale = TRUE ) ) )
}


func_center <- function(
### Mean centers data
### Transforms within columns
df_data
### Data to be transformed
){
  return( data.frame( scale( df_data, center = TRUE, scale = FALSE ) ) )
}


func_asin_sqrt <- function(
### Asin sqrt transform for proportional data exhibiting heteroscedasticity
### Expects propotional data
### Transforms within columns
df_data
### Data to be transformed
){
  return( asin( sqrt( df_data )  ) )
}


######################################
# Filters
######################################

func_filter_downsample <- function(
### Performs downsampling on a full data frame
df_data,
### Data to be down sampled
str_output_dir = NULL,
### Output directory to write results to
...
){
  loginfo( "Downsample counts" )

  # If plotting is needed, make pdf
  f_plotting = !is.null( str_output_dir )
  if( f_plotting )
  {
    pdf( file.path( str_output_dir, "downsampled_populations.pdf" ), useDingbats = FALSE )
  }

  for( i_index in 1:ncol(df_data))
  {
    # Plot before distribution
    if( f_plotting )
    {
      vctr_cur_data = df_data[[ i_index ]]
      barplot( vctr_cur_data, main=paste( "Column", i_index, "(Before and After)"), col = STR_COLOR_BEFORE )
    }
    # Downsample each column (sample)
    df_data[[ i_index ]] = func_downsample( df_data[[ i_index ]] )
    # Plot after distribution
    if( f_plotting )
    {
      barplot( df_data[[ i_index ]], main="", col = STR_COLOR_AFTER, add=TRUE )
      legend( "topright", c("Before","After"), bg = "white", fill = c( STR_COLOR_BEFORE, STR_COLOR_AFTER ))
    }
  }

  if( f_plotting )
  {
    dev.off()
  }
  return( df_data )
}


func_downsample <- function(
### Down samples counts while attempting to keep the distribution shape
### This will not up sample a sample size; in that scenario it will return the original distribution.
### This works on one sample distribution
vctr_i_sample_counts,
### Original feature ( distrbution ) to be downsampled
i_target_population_sum = 100,
### Target total counts in distribution after sampling
...
){
  loginfo( paste( "Downsample to ", i_target_population_sum, " total counts." ) )

  # Get the sum of the population vector
  i_population_sum = sum( vctr_i_sample_counts )
  i_population_length = length( vctr_i_sample_counts )

  # If more counts are requested than exist just return the original distribution
  if( i_target_population_sum > i_population_sum )
  {
    return( vctr_i_sample_counts )
  }

  # 

  # Check and make sure the length of the sample vector is more than
  # one measurement or the sample method will act up
  # If there is only one value in the distribution return it or the amount of sampling requeted
  if( i_population_length == 1 )
  {
    return( min( i_target_population, vctr_i_sample_counts ) )
  }

  # Sample each instance and count instances (positionally)
  vctr_downsampling = rep( 0, i_population_length)
  vi_index_counts = c()
  
  for(i_index in 1:length( vctr_i_sample_counts ))
  {
    vi_index_counts = c(vi_index_counts,rep( i_index, floor(vctr_i_sample_counts[ i_index ] ) ) )
  }
  vi_counts = table( sample( vi_index_counts, size = i_target_population_sum, replace = FALSE ) )
  vctr_downsampling[ as.integer( names( vi_counts ) ) ] = vi_counts
  return( vctr_downsampling )
}


func_filter_by_SD <- function(
### Filter by standard deviation
df_frame,
### Data to filter
str_output_dir = NULL,
### Directory to save plots and such
...
){
  loginfo( "Filter by SD" )  

  # Filter rows by variance
  d_min_sd = 0.1

  # Filter
  vctr_sd = apply( df_frame, 1, function( x ){ sd( x, na.rm = TRUE ) } )
  vctr_sd_no_zero = apply( df_frame, 1, function( x ) { sd( x[ x > 0 ], na.rm = TRUE) } )
  vctr_f_keep = vctr_sd >= d_min_sd

  if( !is.null( str_output_dir ) )
  {
    vctr_colors = rep( STR_COLOR_BEFORE, length( vctr_sd ) )
    vctr_colors[ vctr_f_keep ] = STR_COLOR_STANDARD

    pdf( file.path( str_output_dir, "sd_filter.pdf" ), useDingbats=FALSE )
    plot( sort( vctr_sd ), main = "SD ( with zeros )", col = vctr_colors[ order( vctr_sd ) ] )
    plot( sort( vctr_sd_no_zero ), main = "SD ( with out zeros )", col = vctr_colors[ order( vctr_sd ) ] )
    dev.off()
  }

  df_frame = df_frame[ vctr_f_keep, ]
  return( df_frame )
}


func_filter_by_first_rows <- function(
### Filter by giving the first N rows sorted by expression (head)
df_frame,
### Data to filter
d_top_rows = 50,
### Top number of rows ( in order ) to keep
str_output_dir = NULL,
### Directory to save plots
...
){
  loginfo( "Filter by head (n first top rows )" )

  # Head data frame
  vctr_ordered_by_expression = order( apply( df_frame, 1, sum ), decreasing = TRUE )[ 1:d_top_rows ]
  df_ret = df_frame[ vctr_ordered_by_expression , ]

  # Optionally plot
  if( !is.null( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "head_filter.pdf" ) )
    plt_colors = func_monochromatic_palette()
    # Heatmap before
    func_heatmap( df_frame, str_cor = "euclidean", str_linkage = "average", plt_colors = plt_colors, str_title = "Data before head filter", ... )
    # Heatmap after
    func_heatmap( df_ret, str_cor = "euclidean", str_linkage = "average", plt_colors = plt_colors, str_title = "Data after head filter", ... )
    dev.off() 
  }

  # Return abridged data frame
  return( df_ret )
}


func_filter_by_min_value_occurence <- function(
### Filter by min value in min number samples
df_frame,
### Data to filter
d_min_value = 1,
### Minimum value a sample must have
d_min_occurence = NULL,
### Minimum times the value must occur in a feature to be kept
str_output_dir = NULL,
### Directory to save plots
...
){
  loginfo( "Filter by min value occurence" )

  # If the min occurence is not given, make it 10% of the current data frame size (min 1)
  if( is.null( d_min_occurence ) )
  {
    d_min_occurence = max(1, ceiling( ncol( df_frame ) * .1 ))
  }
  # True of false that a feature passes the filter
  vctr_f_keep = apply( df_frame, 1, function( x ){ return( length( which( x >= d_min_value ) ) >= d_min_occurence ) } )

  # Optionally plot
  if( !is.null( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "value_occurence_filter.pdf" ) )
    plt_colors = func_monochromatic_palette()

    # Calculate values
    vctr_min_value_occurence_keep = apply( df_frame[ vctr_f_keep,], 1, function( x ){ return( length( which( x >= d_min_value ) ) ) })
    vctr_min_value_keep = apply( df_frame[ vctr_f_keep,], 1, function( x ){ return( min( x[x >= d_min_value] ) ) })
    vctr_min_value_occurence = apply( df_frame[ ! vctr_f_keep,], 1, function( x ){ return( length( which( x >= d_min_value ) ) ) })
    vctr_min_value = apply( df_frame[ ! vctr_f_keep,], 1, function( x ){ return( min( x[x >= d_min_value] ) ) })
    i_length_keep = length( vctr_min_value_occurence_keep ) 
    i_length_no_keep = length( vctr_min_value_occurence ) 

    # Plot distributions
    if( ( i_length_keep == 0 ) && ( i_length_no_keep > 0 ) )
    {
       hist( vctr_min_value, col = STR_COLOR_BEFORE,
          main = "Filtering by min value occurence ( All Filtered ) ", xlab = "Min value", ylab = "Occurence across samples" )
    } else if( ( i_length_keep > 0 ) && ( i_length_no_keep == 0 ) )
    {
       hist( vctr_min_value_keep, col = STR_COLOR_STANDARD,
          main = "Filtering by min value occurence ( None Filtered ) ", xlab = "Min value", ylab = "Occurence across samples" )
    } else {
      ret_hist1 = hist( vctr_min_value_keep, plot = FALSE )
      ret_hist2 = hist( vctr_min_value, plot = FALSE )
      d_xlim = c( min( vctr_min_value_keep, vctr_min_value ), max( vctr_min_value_keep, vctr_min_value ) )
      d_ylim = c( 0, max( max( ret_hist1$counts ), max( ret_hist2$counts ) ) )
      vctr_d_breaks = hist( c( vctr_min_value_keep, vctr_min_value ), plot = FALSE )$breaks
      hist( vctr_min_value_keep, col = STR_COLOR_STANDARD, breaks = vctr_d_breaks, xlim = d_xlim, ylim = d_ylim,
          main = "Filtering by min value occurence", xlab = "Min value", ylab = "Occurence across samples" )
      hist( vctr_min_value, col = STR_COLOR_BEFORE, add = TRUE )
    }
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), fill = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ) ) 

    # Plot the values of min value vs min occurence, indicating who was kept
    if( i_length_keep > 0 )
    {
      plot( vctr_min_value_keep, vctr_min_value_occurence_keep, col = STR_COLOR_STANDARD, main = "Filtering by min value occurence", xlab = "Min value", ylab = "Occurence across samples", pch = 16 )
      points( vctr_min_value, vctr_min_value_occurence, col = STR_COLOR_BEFORE, pch = 16 )
    } else {
      plot( vctr_min_value, vctr_min_value_occurence, col = STR_COLOR_BEFORE, main = "Filtering by min value occurence", xlab = "Min value", ylab = "Occurence across samples", pch = 16 )
    }
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), fill = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ) ) 

    # Heatmap before
    func_heatmap( df_frame, str_cor = "euclidean", str_linkage = "average", plt_colors = plt_colors, str_title = "Data before occurence filter", ... )
    if( i_length_keep > 0 )
    {
      # Heatmap after
      func_heatmap( df_frame[ vctr_f_keep,], str_cor = "euclidean", str_linkage = "average", plt_colors = plt_colors, str_title = "Data after occurence filter", ... )
    }
    dev.off() 
  }

  # Return kept data frame
  return( df_frame[ vctr_f_keep, ] )
}


func_filter_by_percentile_percentage <- function(
### Filter rows that do not have a percentage of values in a top percentile of each group
df_frame,
### Data to filter
d_min_percentile = 0.9,
### Minimum percentile of a sample a feature must be in
d_min_occurence = NULL,
### Minimum times the feature must be in the given percentile to be kept
str_output_dir = NULL,
### Directory to save plots
...
){
  loginfo( "Filter by percentile" )

  # Set min occurence if not given
  if( is.null( d_min_occurence ) )
  {
    d_min_occurence = ncol( df_frame ) * 0.1
  }

  # Get the min value for needed percentile for each sample
  vctr_min_values_at_percentile = apply( df_frame, 2, function(x){ return( quantile( x, d_min_percentile ) ) } )

  # Find features that have the min occurence above the given percentile value.
  vctr_f_keep = apply( df_frame, 1, func_min_occurence_at_min_value, vctr_min_values = vctr_min_values_at_percentile ) >= d_min_occurence

  # Optionally plot
  if( !is.null( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "value_percentile_filter.pdf" ) )
    plt_colors = func_monochromatic_palette()

    # Calculate values
    vctr_at_percentile = apply( df_frame, 1, func_min_occurence_at_min_value, vctr_min_values = vctr_min_values_at_percentile )
    vctr_keep = vctr_at_percentile[ vctr_f_keep ]
    vctr_no_keep = vctr_at_percentile[ ! vctr_f_keep ]
    i_length_keep = length( vctr_keep )
    i_length_no_keep = length( vctr_no_keep )

    # Plot distributions
    if( ( i_length_keep == 0 ) && ( i_length_no_keep > 0 ) )
    {
       hist( vctr_at_percentile, col = STR_COLOR_BEFORE,
          main = "Filtering by percentile ( All Filtered ) ", xlab = "Occurence in top percentile across samples", ylab = "Density" )
    } else if( ( i_length_keep > 0 ) && ( i_length_no_keep == 0 ) )
    {
       hist( vctr_at_percentile, col = STR_COLOR_STANDARD,
          main = "Filtering by percentile ( None Filtered ) ", xlab = "Occurence in top percentile across samples", ylab = "Density" )
    } else {
      ret_hist1 = hist( vctr_keep, plot = FALSE )
      ret_hist2 = hist( vctr_no_keep, plot = FALSE )
      d_xlim = c( min( vctr_no_keep, vctr_keep ), max( vctr_keep, vctr_no_keep ) )
      d_ylim = c( 0, max( max( ret_hist1$counts ), max( ret_hist2$counts ) ) )
      vctr_d_breaks = hist( c( vctr_keep, vctr_no_keep ), plot = FALSE )$breaks
      hist( vctr_keep, col = STR_COLOR_STANDARD, breaks = vctr_d_breaks, xlim = d_xlim, ylim = d_ylim,
          main = "Filtering by percentile", xlab = "Occurence in top percentile across samples", ylab = "Density" )
      hist( vctr_no_keep, col = STR_COLOR_BEFORE, add = TRUE )
    }
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), fill = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ) )

    # Heatmap before
    func_heatmap( df_frame, str_cor = "euclidean", str_linkage = "average", plt_colors = plt_colors, str_title = "Data before percentile filter", ... )
    if( i_length_keep > 0 )
    {
      # Heatmap after
      func_heatmap( df_frame[ vctr_f_keep,], str_cor = "euclidean", str_linkage = "average", plt_colors = plt_colors, str_title = "Data after percentile filter", ... )
    }
    dev.off()
  }

  # Return kept data frame
  return( df_frame[ vctr_f_keep, ] )
}

func_min_occurence_at_min_value <- function(
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


func_filter_by_sparsity <- function(
### Filter by percent zeros ( sparsity )
df_frame,
### Data to filter
d_percentile = 0.9,
### Samples must have this max percent of zeros
str_output_dir = NULL,
### Directory to save plots and such
...
){
  loginfo( "Filter by Sparsity" )

  i_nrow = nrow( df_frame )
  i_ncol = ncol( df_frame )

  # Filter rows
  vctr_percent_zero_rows = apply( df_frame, 1, function(x) length( which(x == 0 ) ) ) / i_ncol
  vctr_row_keep = which( vctr_percent_zero_rows < d_percentile )
  vctr_row_col = rep( STR_COLOR_BEFORE, i_nrow )
  vctr_row_col[ vctr_row_keep ] = STR_COLOR_STANDARD

  # Reduce data frame rows
  df_frame = df_frame[ vctr_row_keep, ]

  # Filter columns
  vctr_percent_zero_cols = apply( df_frame, 2, function(x) length( which( x==0 ) ) ) / i_nrow
  vctr_col_keep = which( vctr_percent_zero_cols < d_percentile )
  vctr_col_col = rep( STR_COLOR_BEFORE, i_ncol )
  vctr_col_col[ vctr_col_keep ] = STR_COLOR_STANDARD

  # Filter data frame columns
  df_frame = df_frame[ , vctr_col_keep ]

  # Plot
  if( !is.null( str_output_dir ) )
  {
    pdf( file.path( str_output_dir, "sparsity_filter.pdf" ), useDingbats=FALSE )
    # Plot filtered second so that it is not hidden
    vctr_cur = vctr_percent_zero_rows[vctr_row_keep]
    plot( rnorm( length( vctr_cur ), 1, .1 ), vctr_cur, axes = FALSE, col = STR_COLOR_STANDARD, ylab = "", xlab = "", xlim = c(0,3), pch = 16, cex = .85, ylim = c( min( vctr_percent_zero_rows ) * 0.9, max( vctr_percent_zero_rows ) * 1.1  ) )
    vctr_cur = vctr_percent_zero_rows[ -1 * vctr_row_keep ]
    points( rnorm( length( vctr_cur ), 1, .1 ), vctr_cur, col = STR_COLOR_BEFORE, ylab = "", xlab = "", pch = 16, cex = .85 )

    vctr_cur = vctr_percent_zero_cols[ vctr_col_keep ]
    points( rnorm( length( vctr_cur ), 2, .1 ), vctr_cur, col = vctr_col_col, ylab = "", xlab = "", pch = 16, cex = .85 )
    vctr_cur = vctr_percent_zero_cols[ -1 * vctr_col_keep ]
    points( rnorm( length( vctr_cur ), 2, .1 ), vctr_cur, col = STR_COLOR_BEFORE, ylab = "", xlab = "", pch = 16, cex = .85 )
    boxplot( vctr_percent_zero_rows, vctr_percent_zero_cols, names = c("Row","Columns"), main = "Percent Zeros (columns after row filtering)", ylab = "Percent Zero", add = TRUE )
    legend( "topright", title = "Filtering", legend = c("Filtered","Kept"), col = c( STR_COLOR_BEFORE, STR_COLOR_STANDARD ), pch = 16 ) 
    dev.off()
  }

  return( df_frame ) 
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


func_heatmap <- function(
### Create a heatmap of data frame
df_frame,
### Data for which to show a heatmap
str_cor = "spearman",
### The correlation metric to use
str_linkage = "average",
### Any valid linkage method for hclust
vctr_grouping = NULL,
### Optional Groups the samples
str_output_dir = NULL,
### The directory to output figure and such
plt_colors = NULL,
### Color pallette, if not given ploychromatic default will be used
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


func_cor_matrix <- function(
### Show the correlation of the date in the data frame (samples or features)
df_frame,
### Data for which to show correlation
str_cor = "euclidean",
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
  loginfo( "Visualize correlation" )
  if( f_correlate_samples )
  {
    loginfo( "Sample Correlation" )
  } else {
    loginfo( "Feature Correlation" )
  }

  # Get a consistent color scheme for the plotting
  plt_colors = func_polychromatic_palette( )

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


func_PCA = function(
### Perform ordination with PCA
df_frame,
### Data to ordinate
df_metadata = NULL,
### Metadata to paint the data with
### Features / rows should be in the same order as the data
str_output_dir = NULL,
### Output directory
vctr_str_color_groups = NULL,
...
){
  loginfo( "Ordination by PCA" )

  # Row center and log
  mtrx_scale = t( scale( t( as.matrix( func_log( df_frame ) ) ), center=TRUE, scale=TRUE ) )
  # Remove constant rows
  mtrx_scale = mtrx_scale[ !is.na( mtrx_scale[, colnames( mtrx_scale )[ 1 ] ] ), ]
  # Perfrom PCA
  results_pca = prcomp( mtrx_scale, retx = TRUE )
  # Get percent variance
  vctr_i_variance = results_pca$sdev^2
  vctr_i_percent_variance = ( vctr_i_variance ) / sum( vctr_i_variance )
  i_eigen_count = length( results_pca$sdev )
  # Get the stick break distribution
  result_stick = func_select_pca_components_stick_breaking(vctr_i_percent_variance = vctr_i_percent_variance, str_output_dir = str_output_dir )
  vctr_break_distribution = result_stick$distribution
  i_max_pc = result_stick$count

  # Write loadings
  func_write_frame( results_pca$rotation, file.path( str_output_dir, "pca_ordination_loadings.tsv" ) )
  # write scores
  func_write_frame( results_pca$x, file.path( str_output_dir, "pca_ordination_scores.tsv" ) )

  # Plot in cur directory if output directory not indicated
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }

  # create colorings
  vctr_str_colors = rep( STR_COLOR_STANDARD, ncol(df_frame) )

  # Plot scree vs stick break
  # Plot selected components
  pdf( file.path( str_output_dir, "pca.pdf" ), useDingbats=FALSE )
  for( i_index in 1:i_max_pc )
  {
    plot( results_pca$rotation[, i_index], results_pca$rotation[, i_index + 1 ], 
        xlab = paste("PC", i_index, "(", round( vctr_i_percent_variance[ i_index ]*100, digits = 2 ), "% variance )" ),
        ylab = paste( "PC", i_index + 1, "(", round(vctr_i_percent_variance[ i_index + 1 ]*100, digits = 2 ), "% variance )" ),
        col=vctr_str_colors, main = "Ordination by PCA", pch = 16 )
  }
  vctr_str_comp_colors = rep( "grey",i_eigen_count )
  vctr_str_comp_colors[ 1: i_max_pc ] = "red"
  plot( results_pca, pch = "16", main = "Scree plot", col = vctr_str_comp_colors, xlab = "Component")#, ylab = "Variance", col = vctr_str_comp_colors )

  # Paint with metadata
  for( str_name in names( df_metadata ) )
  {
    vctr_metadata = df_metadata[[ str_name ]]

    vctr_meta_colors = NULL
    vctr_meta_levels = NULL
    vctr_str_colors = NULL

    if( is.numeric( vctr_metadata ) )
    {
      vctr_str_colors = func_polychromatic_palette( length( vctr_metadata ) )
      vctr_meta_levels = c( paste( max(vctr_metadata), "(Max)", sep = " "), paste( min(vctr_metadata), "(Min)", sep = " ") )
      vctr_meta_colors = c( vctr_str_colors[ length( vctr_str_colors ) ], vctr_str_colors[ 1 ] )
    } else {
      if( ! is.factor( vctr_metadata ) )
      {
        vctr_metadata = as.factor( vctr_metadata )
      }
      ret_color_plt = func_factor_to_metadata_color( df_metadata[[ str_name ]] )
      vctr_meta_colors = ret_color_plt$vctr_colors
      vctr_meta_levels = ret_color_plt$vctr_levels
      vctr_str_colors = ret_color_plt$vctr_grouping_colors
    }

    for( i_index in 1:i_max_pc )
    {
      plot( results_pca$rotation[, i_index], results_pca$rotation[, i_index + 1 ], 
        xlab = paste("PC", i_index, "(", round( vctr_i_percent_variance[ i_index ]*100, digits = 2 ), "% variance )" ),
        ylab = paste( "PC", i_index + 1, "(", round(vctr_i_percent_variance[ i_index + 1 ]*100, digits = 2 ), "% variance )" ),
        col=vctr_str_colors, main = paste( "Ordination by PCA (showing metadata ",str_name,")", sep ="" ), pch = 16)
      legend( "topright", vctr_meta_levels, fill = vctr_meta_colors, title = str_name )
    }
  }

  dev.off()
  return( results_pca )
}


func_PCoA <- function(
df_frame,
### Data to ordinate
i_dimensions = 2,
### The number of dimensions to estimate
df_metadata = NULL,
### Data rame of metadata with features ordered as the data data frame
### If given, ordinations will be painted with these metadata
str_metric = "bray",
### The distance metric to use (must be given or mtrx_dist must be given)
### "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard",
### "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
mtrx_dist = NULL,
### The distance matrix to use (must be given or str_metric must be given)
str_output_dir = NULL,
### Output directory
...
){
  ###
  # Note:
  # Much thanks to the following post which inspired much of this code. Well done!
  # http://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
  ##

  loginfo( paste( "NMDS (PCoA) Ordination using", str_metric, "distance/dissimilarity", sep = " " ) )

  # Plot ordination
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }

  # Start pdf
  pdf( file.path( str_output_dir, paste( "pcoa_", str_metric,".pdf", sep = "" ) ), useDingbats=FALSE )

  # Check that we can make or have a distance matrix
  if( is.null( mtrx_dist ) )
  {
    if( is.null( str_metric ) )
    {
      logerror( "NMDS (PCoA) Ordination. Please give the either a distance/dissimilarity metric or matrix. Ordination was not performed." )
      return( FALSE )
    }

    # Generate a distance matrix
    mtrx_dist = func_get_custom_distance_matrix( t( df_frame ), str_metric )
  }

  # Perform NMDS
  ret_NMDS = metaMDS( mtrx_dist, k=i_dimensions, autotransfor = FALSE, trymax = 100 )

  # Plot ordination but only samples
  # Plot sample ordination with sample names
  for( str_plot_type in c( "points", "text" ) )
  {
    ordiplot( ret_NMDS, choices = c(1,2), display = "sites", type = str_plot_type, xlab = "NMDS 1", ylab = "NMDS 2",
            main = paste( "NMDS Ordination (",str_metric,") Stress=", ret_NMDS$stress, sep = " " ))
  }

  # Plot stress
  stressplot( ret_NMDS )


  # If no metadata is given return at this point
  if( is.null( df_metadata ) )
  {
    dev.off()
    return( TRUE )
  }

  # Plot sample ordination painted by metadata
  for( str_name in names( df_metadata ) )
  {
    vctr_metadata = df_metadata[[ str_name ]]
    if( is.numeric( vctr_metadata ) )
    {
      # If continuous, plot with gradient
      ordisurf( ret_NMDS, vctr_metadata, choices = c(1,2), col = STR_COLOR_AFTER, xlab = "NMDS1", ylab = "NMDS2",
                main = paste( "NMDS Ordination(",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      plot( ret_NMDS$points[,1], ret_NMDS$points[,2], col = STR_COLOR_STANDARD, add = TRUE )
      legend( "topright", c( paste( max(vctr_metadata), "(Max)", sep = " " ), paste( min(vctr_metadata), "(Min)", sep = " " )),
              title = str_name ) 

      # Plot samples as names
      ordisurf( ret_NMDS, vctr_metadata, choices = c(1,2), col = STR_COLOR_AFTER, xlab = "NMDS1", ylab = "NMDS2",
                main = paste( "NMDS Ordination(",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      orditorp( ret_NMDS, choices = c(1,2), display = "sites", col = STR_COLOR_STANDARD )
      legend( "topright", c( paste( max(vctr_metadata), "(Max)", sep = " " ), paste( min(vctr_metadata), "(Min)", sep = " " )),
              title = str_name ) 

    } else {
      # If not numeric make factor and plot
      if( !is.factor( vctr_metadata ) )
      {
        vctr_metadata = as.factor( vctr_metadata )
      }

      # Make color
      ret_meta_color = func_factor_to_metadata_color( vctr_metadata )
      vctr_metadata_levels = ret_meta_color$vctr_levels
      vctr_plt_colors = ret_meta_color$vctr_colors
      vctr_colors = ret_meta_color$vctr_grouping_colors

      # If discrete data, plot with ellipse
      plot( ret_NMDS$points[,1], ret_NMDS$points[,2], col = vctr_colors, add = TRUE,
            main = paste( "NMDS Ordination (",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )

      # If discrete data, plot with ellipse
      plot( ret_NMDS$points[,1], ret_NMDS$points[,2], col = vctr_colors, add = TRUE,
            main = paste( "NMDS Ordination (",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )

      # Plot each group to change color for groups
      for( i_index in 1:length( vctr_colors ) )
      {
        str_cur_level = vctr_metadata_levels[ i_index ]
        str_cur_color = vctr_plt_colors[ i_index ]
        ordiellipse( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", kind = "se", conf = 0.95 ,
                   draw = "line", w = NULL, col = str_cur_color, label = FALSE, show.groups = str_cur_level )
        ordispider( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", col = str_cur_color, label = FALSE, show.groups = str_cur_level )
      }
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )

      # Plot samples as names
      ordiplot( ret_NMDS, choices = c(1,2), type = "none", xlab = "NMDS1", ylab = "NMDS2",
                main = paste( "NMDS Ordination (",str_metric,") Stress=",ret_NMDS$stress, sep = " "  ) )

      # Plot each group to change color for groups
      for( i_index in 1:length( vctr_colors ) )
      {
        str_cur_level = vctr_metadata_levels[ i_index ]
        str_cur_color = vctr_plt_colors[ i_index ]
        ordiellipse( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", kind = "se", conf = 0.95 ,
                   draw = "line", w = NULL, col = str_cur_color, label = FALSE, show.groups = str_cur_level )
        ordispider( ret_NMDS, vctr_metadata, choices = c(1,2), display = "sites", col = str_cur_color, label = FALSE, show.groups = str_cur_level )
      }
      orditorp( ret_NMDS, choices = c(1,2), display = "sites", col = vctr_colors )
      legend( "topright", vctr_metadata_levels, fill = vctr_plt_colors, title = str_name )
    }
  }

  dev.off()
  return( TRUE )
}


############################################
# Quality Control
############################################

func_QC <- function(
### Perform quality control on matrix of data
df_frame,
### Data to have quality control on
df_metadata = NULL,
### Metadata to use for correction
func_handle_outliers = func_do_nothing,
### Function to handle outliers in metadata
func_imputation = func_do_nothing,
### Function for imputing data
func_batch_correct = func_do_nothing,
### The method for batch correction
str_output_dir = NULL
###
){
  loginfo( "func_QC" )

  df_cleaned_data = df_frame
  df_cleaned_metadata = df_metadata

  # Handle outliers
  # TODO
  
  # Handle imputation
  if( !is.null( df_metadata ) && ( ! is.null( func_imputation ) ) ) 
  {
    list_ret = func_imputation( df_frame, df_metadata )
    df_cleaned_data = list_ret$data
    df_cleaned_metadata = list_ret$metadata
  }

  # Handle batch correction
  # TODO
  return( list( data=df_cleaned_data, metadata=df_cleaned_metadata ) )
}


#--- Metadata imputation


func_complete_case <- function(
### Handle missing data by reducing the metadata and data observations
### to those that do not have missing data
df_data,
### Data to reduce based on metadata (assumed to be the same observations as metadata)
df_metadata
### Data to reduce to observation with no missing values
){
  loginfo( "Handle missing values by reducing metadata and data to observations without missing values" )
  df_complete_metadata = na.omit( df_metadata )
  df_complete_data = df_data[ ,row.names( df_complete_metadata ) ]
  return( list( data = df_complete_data, metadata = df_complete_metadata ) )
}


func_impute_mean <- function(
### Impute missing values with a mean for numerics and integers
### binary and factor data make a sperate NA category
df_data,
### Data measurements (will also be imputed)
df_metadata
### Metadata measurements to impute missing data
){
  loginfo( "Impute with mean" )

  # Reset the metadata
  df_metadata = func_impute_with_center_helper( df_metadata, mean )
  # Reset the data
  df_data = func_impute_with_center_helper( df_data, mean )
  # Return
  return( list( data=df_data, metadata=df_metadata ) )
}


func_impute_median <- function(
### Impute missing values with a median for numerics and integers
### binary and factor data make a sperate NA category
df_data,
### Data measurements (will also be imputed)
df_metadata
### Metadata measurements to impute missing data
){
  loginfo( "Impute with median" )

  # Reset the metadata
  df_metadata = func_impute_with_center_helper( df_metadata, median )
  # Reset the data
  df_data = func_impute_with_center_helper( df_data, median )
  # Return
  return( list( data=df_data, metadata=df_metadata ) )
}


func_impute_with_center_helper <- function(
df_matrix,
### Values to impute
func_center
### Function to use to make the value to imput with
){
  # Reset the data
  for( str_feature in names( df_matrix ) )
  {
    vctr_values = df_matrix[[str_feature]]
    vctr_omit = union( which( is.na( vctr_values ) ),which( is.nan( vctr_values ) ) )
    if( is.numeric( vctr_values ) )
    {
      df_matrix[ vctr_omit, str_feature ] = func_center( vctr_values, na.rm = TRUE )
    } else if ( is.integer( vctr_values ) )
    {
      df_matrix[ vctr_omit, str_feature ] = round( func_center( vctr_values, na.rm = TRUE ) )
    } else if ( is.factor( vctr_values ) || is.ordinal( vctr_values ) )
    {
      levels( vctr_values ) = c( levels(vctr_values), "NA" )
      vctr_values[ vctr_omit ] = "NA"
      df_matrix[, str_feature] = vctr_values
    }
  }
  return( df_matrix )
}


# TODO
# func_impute_multiple <- function(
#


#---- Batch affect correction


func_PCA_skimming <- function(
### Removes batch affect by removing PCs correclated with the metadata
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_PCA_skimming" )
  df_cleaned = df_frame
  return( df_cleaned )
}


func_combat <- function(
### Removes batch affect by using the combat method
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_combat" )

  df_cleaned = df_frame
  return( df_cleaned )
}


func_control_for_batching <- function(
### Removes batch affects by regession
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_control_for_batching" )

  df_cleaned = df_frame
  retrun( df_cleaned )
}


func_edgeR_batch_control <- function(
### Removes batch affects with edgr
df_frame,
### Data to correct
df_metadata
### Metadata to control for
){
  loginfo( "func_edgeR_batch_control" )
  library( edgeR )

  print( removeBatchEffect( x = df_frame, covariates = df_metadata ) )
}


###########################################
# Subset selection
###########################################

func_select_by_ANOVA <- function(
### Select genes by anova or t-test
df_frame,
### Data to perform selection from
vctr_groupings,
### Data to make groupings
str_output_dir
### Output directory
){
  loginfo( "Feature selection by ANOVA" )
  i_group_count = length( unique( vctr_groupings ) )
  list_results = list()
  if( i_group_count == 1 )
  {
    print( paste( "Can not perform or t-test, number of groups found are ", i_group_count, "." ) )
  } else if( i_group_count == 2 ) {
    ### Perform T-test
  } else if ( i_group_count > 2 ) {
    ### Perform ANOVA
    for( row_name in row.names( df_frame ) )
    {
      results_anova = aov( as.vector( unlist(df_frame[ row_name, ])) ~ vctr_groupings )
      list_results[[ row_name ]] = summary(results_anova)[[1]]["vctr_groupings",]
    }
  }
  return( list_results )
}


func_select_features_by_pca <- function(
### Select by PCA
df_frame,
### Data to perform selection from
i_PC_Components = NULL,
### Number of PC components
i_count = NULL,
### The number of top ranked features to select
str_output_dir = NULL
### Output directory to document to
){
  loginfo( "Feature selection by PCA" )
  # Row center and log
  mtrx_scale = t( scale( t( as.matrix( func_log( df_frame ) ) ), center=TRUE, scale=TRUE ) )
  # Remove constant rows
  mtrx_scale = mtrx_scale[ !is.na( mtrx_scale[, colnames( mtrx_scale )[ 1 ] ] ), ]
  # Perfrom PCA
  results_pca = princomp( mtrx_scale, cor = TRUE )
  # For each component get the top genes
  mtrx_scores = results_pca$scores
  # How many genes exist
  i_gene_counts = nrow( mtrx_scores )
  # Get percent variance
  vctr_i_variance = results_pca$sdev^2
  vctr_i_percent_variance = ( vctr_i_variance ) / sum( vctr_i_variance )
  i_eigen_count = length( results_pca$sdev )

  # If check output directory path
  if( is.null( str_output_dir ) )
  {
    str_output_dir = "."
  }

  # If no count is given us the upper quartile expressed
  if( is.null( i_count ) )
  {
    i_count = floor( i_gene_counts * 0.25 )
  }

  # Get the stick break distribution
  if( is.null( i_PC_Components ) )
  {
    result_stick = func_select_pca_components_stick_breaking( vctr_i_percent_variance = vctr_i_percent_variance, str_output_dir = str_output_dir )
    i_PC_Components = result_stick$count
  }

  # Holds the selection per component
  vctr_return_feature_list = c()

  # For each component select the top genes
  for( i_index_component in 1:i_PC_Components )
  {
    vctr_genes = row.names( mtrx_scores[ order( abs( mtrx_scores[,i_index_component] ) ), ] )[ 1 : i_count ]
    if( ! is.null( str_output_dir ) )
    {
      print( str_output_dir )
      func_write_frame( vctr_genes, file.path( str_output_dir, paste( "Extreme",i_count, "genes","by","PCA_component",i_index_component,".txt", sep="_" ) ) )
    }
    vctr_return_feature_list = unique( c( vctr_return_feature_list, vctr_genes ) )
  }
  return( vctr_return_feature_list )
  ### Returns a list of selected genes
}


###########################################
# Functional analysis
###########################################

func_david <- function(
### Perform functional analysis by interfacting with the david website
vctr_genes,
### Genes to evaluate
str_output_dir
### Output directory
){
  library("RDAVIDWebService")

  david_service = DAVIDWebService$new(email=str_email)
  david_list = addList( david_service, vctr_genes, idType=str_annotation_type, listName= paste( format(Sys.time(), "%b_%d_%Y_%X"), "automated_analysis",sep="_"), listType="Gene")
  setAnnotationCategories( david_service, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
  termCluster = getClusterReport( david_service, type="Term" )
  getClusterReportFile( david_service, type="Term", fileName=path.join( str_output_dir, "termClusterReport1.tab" ) )
  
}


func_gsea_preranked <- function(
### Performs GSEA with preranked gene lists

){
# From the limma package
wilcoxGST()
}


func_gsea_response <- function(
### Performs GSEA with response

){
GSEABase
}


func_stringdb <- function(
### Interface with String DB to make gene networks
vctr_genes,
### Genes to make a network with
str_output_dir
### Output directory
){

}


###########################################
# Generic visualizations
###########################################


#func_matrix_plot <- function(
#### This function plots raw dataframe
#df_frame,
#### Data frame oto plot
#title = ""
#### Title to figure
#){
#  heatmap.3b( df_frame, dendrogram = "both", scale = "center", col = func_monochromatic_palette( ), 
#              trace = "none", tracecol = C_STR_DENSITY_COLOR, density.info = "histogram", main = title )
#}


###########################################
# Module Infrastructure
###########################################

func_do_nothing <- function(
### This function does nothing
df_frame,
### Data frame to act on
df_metadata = null
### Metadata data frame
){
  return( df_frame )
  ### Returns the identical data frame
}


func_get_batch_module <- function(
### This function take a keyword and returns the function
### that would perform batch analysis associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    combat = func_combat,
    edgeR = func_edgeR_batch_control,
    none, func_do_nothing,
    pca_skimming = func_PCA_skimming,
    regression = func_control_by_regression )
  # Return function
  return( func_ret )
  ### Returns batch function to perform the needed analysis 
}


func_get_discriminant_module <- function(
### This function take a keyword and returns the function
### that would perform that discriminate feature analysis associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    anova = func_discriminate_by_anova,
    none = func_do_nothing )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed discriminant feature analysis 
}


func_get_filter_module <- function(
### This function take a keyword and returns the function
### that would perform the filtering associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    downsample = func_filter_downsample,
    head = func_filter_by_first_rows,
    none = func_do_nothing,
    occurence = func_filter_by_min_value_occurence,
    percentile = func_filter_by_percentile_percentage,
    sd = func_filter_by_SD,
    sparsity = func_filter_by_sparsity )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed filtering 
}


func_get_impute_module <- function(
### This function take a keyword and returns the function
### that would perform the imputation associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    complete_case = func_complete_case,
    mean = func_impute_mean,
    median = func_impute_median,
    multiple_imputation = func_impute_multiple,
    none = func_do_nothing )
  # Return function
  return( func_ret )
  ### Returns function to perform the imputation needed
}


func_get_normalization_module <- function(
### This function taks a keyword and returns the function
### that woudl perform the imputation associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    cpm = func_cpm,
    none = func_do_nothing,
    upper = func_upper_quartile )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed ordination
}


func_get_ordinate_module <- function(
### This function take a keyword and returns the function
### that would perform the ordination associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    none = func_do_nothing,
    pca = func_PCA,
    pcoa = func_PCoA )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed ordination
}


func_get_select_sample_groups_module <- function(
### This function take a keyword and returns the function
### that would perform the sample group selection associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    anova = func_select_by_anova,
    none = func_do_nothing )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed sample group selection 
}


func_get_select_feature_groups_module <- function(
### This function take a keyword and returns the function
### that would perform the feature group associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    none = func_do_nothing,
    pca = func_select_features_by_pca )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed feature group selection 
}


func_get_transform_module <- function(
### This function take a keyword and returns the function
### that would perform the data transform associated with the keyword
str_key
### Key word to act on
){
  # Normalize keyword
  if ( is.null( str_key ) )
  {
    return( NULL )
  }
  str_key = tolower( str_key )
  # Get function
  func_ret = switch( str_key,
    asin_sqrt = func_asin_sqrt,
    center = func_center,
    log = func_log,
    none = func_do_nothing,
    tss = func_tss,
    zscale = func_zscale )
  # Return function
  return( func_ret )
  ### Returns function to perform the needed data transform
}
