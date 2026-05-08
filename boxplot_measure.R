boxplot_measure <- function(nsimul,measure,name,yrange=c(0,1)){
  
  library(tidyverse)
  library(hrbrthemes)

  data <- data.frame(method=c(rep("only U",nsimul),
                              rep("only Z",nsimul),
                              rep("early fusion",nsimul),
                              rep("rho=0.25",nsimul),
                              rep("rho=0.5",nsimul),
                              rep("rho=0.75",nsimul)),
                     measure=c(measure$`only U`,
                               measure$`only Z`,
                               measure$`early fusion`,
                               measure$`rho=0.25`,
                               measure$`rho=0.5`,
                               measure$`rho=0.75`)
  )
  
  data$method <- factor(data$method, 
                 levels=c("only U", "only Z", "early fusion", 
                          "rho=0.25", "rho=0.5", "rho=0.75"))
  
  my_palette <- c(
    "only U" = "#4477AA",     
    "only Z" = "#66CCEE",     
    "early fusion" = "#228833",
    "rho=0.25" = "#CCBB44",    
    "rho=0.5" = "#AA3377",   
    "rho=0.75" = "#EE6677"     
  )
  
  # # Plot
   g <- data %>%
     ggplot(aes(x=method, y=measure, fill=method)) +
     geom_boxplot() +
     scale_fill_manual(values=my_palette) +
     theme_ipsum() +
     theme(
       legend.position="none",
       plot.title = element_text(size=11),
       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)
     ) +
     ggtitle(name) +
     xlab("") +
     ylab("")
  
  
  gg <- g + 
    scale_x_discrete(labels = c("only U","only Z","early fusion",
                                expression(rho ~ "= 0.25"),
                                expression(rho ~ "= 0.50"),
                                expression(rho ~ "= 0.75")))+
    scale_y_continuous(limits=yrange)
    

}

