library(ggplot2)
library(reshape2)

plot_heatmap <- function(input_matrix, legend_title, has.legend=T, 
                         colormap="viridis", min=0, max=3){
  
  rownames(input_matrix) <- NULL
  colnames(input_matrix) <- NULL
  
  long_matrix <- melt(input_matrix)
  base_plot <- ggplot(data = long_matrix, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() + labs(fill=legend_title) + 
    scale_fill_viridis_c(option = colormap, limits=c(min, max))+
    theme_minimal()
  
  if (has.legend){
    final_plot <- base_plot 
  } else{
    final_plot <- base_plot + theme(legend.position = "none")
  }
  
  return(final_plot)
}
