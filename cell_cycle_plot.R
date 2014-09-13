#######################
# Get default ids for cell cycle markers
#######################

func_get_default_G0_Ids <- function(
### Get default ids
){
  return( c("G0_1") )
}

func_get_default_G1_Ids <- function(
### Get default ids
){
  return( c("G1_1","G1_2") )
}

func_get_default_S_Ids <- function(
### Get default ids
){
  return( c("S_1","S_2","S_3") )
}

func_get_default_G2_Ids <- function(
### Get default ids
){
  return( c("G2_1","G1_2") )
}

func_get_default_prophase_Ids <- function(
### Get default ids
){
  return( c("Prophase_1") )
}

func_get_default_prometaphase_Ids <- function(
### Get default ids
){
  return( c("Prometaphase_1","Prometaphase_2") )
}

func_get_default_metaphase_Ids <- function(
### Get default ids
){
  return( c("Metaphase_1","Metaphase_2","Metaphase_3") )
}

func_get_default_anaphase_Ids <- function(
### Get default ids
){
  return( c("Anaphase_1","Anaphase_2") )
}

func_get_default_telophase_Ids <- function(
### Get default ids
){
  return( c("Telophase_1") )
}

#####################
# Utility methods
####################

func_circle_xy <- function(
### Returns the x/y coordiantes for the end point of a line from origin
### with an angle and radius
d_angle,
### Angle
d_radius
### Radius
){
  return( c( d_radius * cos( ( d_angle * pi )/180 ), d_radius * sin( ( d_angle * pi )/180 )))
}


func_plot_circle_segment <- function(
### Plots a spider graph segment
d_bar_width,
### How wide the angle is ( how many degrees )
i_border_points,
### The number points to use to make segments; the more the smoother the plot
d_angle,
### The angle of the segment ( where the segment starts on the circle )
d_radius,
### The radius of the segment ( how long it is )
f_solid_segment = TRUE,
### True indicates the segment will be plotted solid, otherwise plotted as a line
str_color = "Black",
str_border = "Black"
){

  # Initialize
  vctr_points_x = c(0)
  vctr_points_y = c(0)
  d_smoothing_angle = d_angle - ( d_bar_width / 2.0 )
  d_increment = d_bar_width / i_border_points

  # Create points, the more the smoother
  # This is controlled by d_border_ponts
  for( i_point in 1:i_border_points )
  {
    vi_points = func_circle_xy( d_smoothing_angle, d_radius )
    vctr_points_x = c( vctr_points_x, vi_points[ 1 ] )
    vctr_points_y = c( vctr_points_y, vi_points[ 2 ] )
    d_smoothing_angle = d_smoothing_angle + d_increment
  }

  # end with origin
  vctr_points_x = c( vctr_points_x, 0 )
  vctr_points_y = c( vctr_points_y, 0 )

  if( f_solid_segment )
  {
    polygon( vctr_points_x, vctr_points_y, col = str_color, border = str_border )
  } else {
    lines( vctr_points_x, vctr_points_y, col = str_color )
  }
}


func_get_group_expression <- function(
df_data,
vctr_group
){
  if( length( vctr_group ) > 1 )
  {
    return( median( apply( df_data[ vctr_group,  ], 1, function( x ){ mean( x, trim = 0.05 ) } ) ))
  }
  return ( mean( as.matrix( df_data[ vctr_group,  ], trim = 0.05 ) ) )
}


##########################
# Main script
##########################


func_cell_cycle_plot <- function(
### Plots the trimmed mean expression of transcripts in different cell cycle phases
df_data,
### The data to be plotted ( rows = transcripts, columns = samples )
vctr_str_G0 = NULL,
### Row ids for G0 phase associated transcripts
vctr_str_G1 = NULL,
### Row ids for G1 phase associated transcripts
vctr_str_S = NULL,
### Row ids for S phase associated transcripts
vctr_str_G2 = NULL,
### Row ids for G2 phase associated transcripts
vctr_str_prophase = NULL,
### Row ids for Prophase phase associated transcripts
vctr_str_prometaphase = NULL,
### Row ids for Prometaphase phase associated transcripts
vctr_str_metaphase = NULL,
### Row ids for Metaphase phase associated transcripts
vctr_str_anaphase = NULL,
### Row ids for Telophase phase associated transcripts
vctr_str_telophase = NULL
){
  # Colors
  plt_cur = func_polychromatic_palette()
  i_color_count = length( plt_cur )

  # Get fold change by mean
  df_data_centered = df_data / mean( as.matrix( df_data ), trim = 0.05 )

  # Get the min and max number in the data for scaling
  # One the off chance that everything is the same number
  # Set min and max so that the numbers end up giving a chart
  # Indicating the values are from the center of the data
  d_min = min( df_data_centered )
  d_max = max( df_data_centered )
  if( d_max == d_min )
  {
    d_max = d_max * 2
    d_min = 0
  }

  # Get the feature ids if not given
  if( is.null(vctr_str_G0) ){ vctr_str_G0 = func_get_default_G0_Ids() }
  if( is.null(vctr_str_G1) ){ vctr_str_G1 = func_get_default_G1_Ids() }
  if( is.null(vctr_str_S) ){ vctr_str_S = func_get_default_S_Ids() }
  if( is.null(vctr_str_G2) ){ vctr_str_G2 = func_get_default_G2_Ids() }
  if( is.null(vctr_str_prophase) ){ vctr_str_prophase = func_get_default_prophase_Ids() }
  if( is.null(vctr_str_prometaphase) ){ vctr_str_prometaphase = func_get_default_prometaphase_Ids() }
  if( is.null(vctr_str_metaphase) ){ vctr_str_metaphase = func_get_default_metaphase_Ids() }
  if( is.null(vctr_str_anaphase) ){ vctr_str_anaphase = func_get_default_anaphase_Ids() }
  if( is.null(vctr_str_telophase) ){ vctr_str_telophase = func_get_default_telophase_Ids() }

  # Get expression data for groups
  d_G0 = func_get_group_expression( df_data_centered, vctr_str_G0 )
  d_G0 = ( d_G0 - d_min ) / ( d_max - d_min )
  d_G1 = func_get_group_expression( df_data_centered, vctr_str_G1 )
  d_G1 = ( d_G1 - d_min ) / ( d_max - d_min )
  d_S = func_get_group_expression( df_data_centered, vctr_str_S )
  d_S = ( d_S - d_min ) / ( d_max - d_min )
  d_G2 = func_get_group_expression( df_data_centered, vctr_str_G2 )
  d_G2 = ( d_G2 - d_min ) / ( d_max - d_min )
  d_prophase = func_get_group_expression( df_data_centered, vctr_str_prophase )
  d_prophase = ( d_prophase - d_min ) / ( d_max - d_min )
  d_prometaphase = func_get_group_expression( df_data_centered, vctr_str_prometaphase )
  d_prometaphase = ( d_prometaphase - d_min ) / ( d_max - d_min )
  d_metaphase = func_get_group_expression( df_data_centered, vctr_str_metaphase )
  d_metaphase = ( d_metaphase - d_min ) / ( d_max - d_min )
  d_anaphase = func_get_group_expression( df_data_centered, vctr_str_anaphase )
  d_anaphase = ( d_anaphase - d_min ) / ( d_max - d_min )
  d_telophase = func_get_group_expression( df_data_centered, vctr_str_telophase )
  d_telophase = ( d_telophase - d_min ) / ( d_max - d_min )

  # Text constants
  str_G0 = "G0"
  str_G1 = "G1"
  str_S = "S"
  str_G2 = "G2"
  str_prophase = "Prophase"
  str_prometaphase = "Prometaphase"
  str_metaphase = "Metaphase"
  str_anaphase = "Anaphase"
  str_telophase = "Telophase"

  # Make new plot
  plot( 0, 0, ylim = c(-12,12), xlim = c(-12,12), main = "Expression by Cell Cycle", yaxt = "n", xaxt = "n", ann = FALSE, frame.plot = FALSE )

  # Plot G0
#  func_plot_circle_segment( d_bar_width = 1, i_border_points = 100, d_radius = 4 + ( 4 * d_G0 ), d_angle = 1, str_color = plt_cur[ floor( i_color_count * d_G0 ) ] )
#  text( str_G0 )

  # Plot G1
  func_plot_circle_segment( d_bar_width = 55, i_border_points = 100, d_radius = 4 + ( 4 * d_G1 ), d_angle = 300, str_color = plt_cur[ floor( i_color_count * d_G1 ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 300, d_radius = 9 )
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_G1 )

  # Plot S
  func_plot_circle_segment( d_bar_width = 138, i_border_points = 100, d_radius = 4 + ( 4 * d_S ), d_angle = 204, str_color = plt_cur[ floor( i_color_count * d_S ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 204, d_radius = 9 )
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_S )


  # Plot G2
  func_plot_circle_segment( d_bar_width = 60, i_border_points = 100, d_radius = 4 + ( 4 * d_G2 ), d_angle = 104, str_color = plt_cur[ floor( i_color_count * d_G2 ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 104, d_radius = 9 )
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_G2 )

  # Plot Prophase
  func_plot_circle_segment( d_bar_width = 20, i_border_points = 100, d_radius = 4 + ( 4 * d_prophase ), d_angle = 63, str_color = plt_cur[ floor( i_color_count * d_prophase ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 63, d_radius = 9 )
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_prophase )


  # Plot Prometaphase
  func_plot_circle_segment( d_bar_width = 20, i_border_points = 100, d_radius = 4 + ( 4 * d_prometaphase ), d_angle = 0, str_color = plt_cur[ floor( i_color_count * d_prometaphase ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 0, d_radius = 10.5 )
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_prometaphase )

  # Plot Metaphase
  func_plot_circle_segment( d_bar_width = 20, i_border_points = 100, d_radius = 4 + ( 4 * d_metaphase ), d_angle = 21, str_color = plt_cur[ floor( i_color_count * d_metaphase ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 21, d_radius = 10 )
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_metaphase )


  # Plot Anaphase
  func_plot_circle_segment( d_bar_width = 20, i_border_points = 100, d_radius = 4 + ( 4 * d_anaphase ), d_angle = 42, str_color = plt_cur[ floor( i_color_count * d_anaphase ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 42, d_radius = 10 )
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_anaphase )


  # Plot Telophase
  func_plot_circle_segment( d_bar_width = 20, i_border_points = 100, d_radius = 4 + ( 4 * d_telophase ), d_angle = 339, str_color = plt_cur[ floor( i_color_count * d_telophase ) ] )
  vctr_d_coordinates = func_circle_xy( d_angle = 339, d_radius = 10 ) 
  text( vctr_d_coordinates[1], vctr_d_coordinates[2], str_telophase )

  # Plot center
  func_plot_circle_segment( d_bar_width = 365, i_border_points = 360, d_radius = 4, d_angle = 0, str_color = "black", str_border = "black" )
  func_plot_circle_segment( d_bar_width = 365, i_border_points = 360, d_radius = 3.9, d_angle = 0, str_color = "white", str_border = "white" )
}

