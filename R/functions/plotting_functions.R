# R/functions/plotting_functions.R

#' Create a volcano plot from limma results
#'
#' @param de_results A data frame from limma's topTable function, must contain
#'                   'logFC', 'adj.P.Val', and 'Gene' columns.
#' @param p_threshold Adjusted p-value threshold for significance.
#' @param fc_threshold Log2 fold change threshold for significance.
#'
#' @return A ggplot object representing the volcano plot.
#'
plot_volcano_r <- function(de_results, p_threshold = 0.05, fc_threshold = 1.0) {
  
  # Add a column to categorize proteins (Upregulated, Downregulated, Not Significant)
  de_results <- de_results %>%
    mutate(significance = case_when(
      adj.P.Val < p_threshold & logFC > fc_threshold  ~ "Upregulated",
      adj.P.Val < p_threshold & logFC < -fc_threshold ~ "Downregulated",
      TRUE                                            ~ "Not Significant"
    ))
  
  # Create the plot
  p <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    theme_bw(base_size = 14) +
    labs(
      title = "Volcano Plot of Differential Protein Abundance",
      x = "Log2 Fold Change",
      y = "-log10(Adjusted P-value)"
    ) +
    # Add threshold lines
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
    # Add labels to the most significant proteins
    geom_text_repel(
      data = head(de_results %>% filter(significance != "Not Significant") %>% arrange(adj.P.Val), 10),
      aes(label = Gene),
      size = 3.5,
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    theme(legend.position = "top")
    
  return(p)
}
