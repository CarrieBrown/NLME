plot_parm_dist <- function(subset, parm, title, xlab, file, true_value){
  ggplot(subset, aes(x=x, color=method, fill=method)) + geom_density(alpha=.2) +
    labs(title = title) +
    ylab("Density") + xlab(xlab) + 
    scale_color_manual(values = c("#249ab5", "#a5228d")) +
    scale_fill_manual(values = c("#249ab5", "#a5228d")) +
    geom_vline(xintercept=true_value) + theme(legend.position = "none",
                                              plot.title = element_text(hjust = 0.5),
                                              plot.subtitle = element_text(hjust = 0.5))
  
  ggsave(paste0(file,"_distribution.png"), width = 6, height = 5)
}

plot_parm_ci <- function(subset, parm, file, true_value){
  iml_subset <- subset %>% filter(method=="iml") %>% arrange(upper+((upper-lower)/2))
  iml_subset$id <- 1:nrow(iml_subset)
  nlm_subset <- subset %>% filter(method=="nlm") %>% arrange(upper+((upper-lower)/2))
  nlm_subset$id <- 1:nrow(nlm_subset)
  
  ggplot(iml_subset, aes(x=id, y=x, ymin=lower, ymax=upper, color=capture)) +
    geom_linerange()+
    coord_flip() + 
    labs(title="Quadrature") +
    geom_hline(yintercept = true_value) +
    scale_color_manual(values = c("#c7c8ca", "#a5228d", "#249ab5")) + 
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(expand=c(0,0))
    
    ggsave(paste0(file, "_iml_ci.png"), width = 4, height = 15)
  
  ggplot(nlm_subset, aes(x=id, y=x, ymin=lower, ymax=upper, color=capture)) +
    geom_linerange()+
    coord_flip() + 
    labs(title="Pseudo-likelihood") +
    geom_hline(yintercept = true_value) +
    scale_color_manual(values = c("#c7c8ca", "#a5228d", "#249ab5")) + 
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(expand=c(0,0))
    
    ggsave(paste0(file,"_nlm_ci.png"), width = 4, height = 15)
}
