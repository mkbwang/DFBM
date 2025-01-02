library(ggplot2)
library(reshape2)

plot_heatmap <- function(input_matrix, legend_title, has.legend=T,
                         ynames="Gene", xnames="Sample",
                         color.low="black", color.high="white", color.na="#6ca0f5",
                         min=0, max=1){

  rownames(input_matrix) <- NULL
  colnames(input_matrix) <- NULL

  long_matrix <- melt(input_matrix)
  base_plot <- ggplot(data = long_matrix, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() + labs(x = xnames, y=ynames, fill=legend_title) +
    scale_fill_gradient(low=color.low, high=color.high, na.value=color.na)+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())

  if (has.legend){
    final_plot <- base_plot
  } else{
    final_plot <- base_plot + theme(legend.position = "none")
  }

  return(final_plot)
}
