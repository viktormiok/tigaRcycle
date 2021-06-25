#########################################################################################
#
#       function for plotting individual gene over time
#
#########################################################################################

circadian_gene_plot <- function(gene_name,
                                data1,
                                tp1,
                                group1,
                                data2,
                                tp2,
                                group2){
                    # gene_name: Name of the gene of interest, 
                    #            which exist in the row names of the data
                    # data1: First data set
                    # tp1: Time points of the first data set
                    # group1: Group information of the first data set
                    # data2: Secon data set
                    # tp1: Time points of the first dat set
                    # group2: Group information of the second data set
                    if (!is(gene_name, "character")) {
                      stop("Input (gene_name) is of wrong class.")
                    }
                    if (!is(data1, "matrix")) {
                      stop("Input (data1) is of wrong class.")
                    }
                    if (!is(tp1, "factor")) {
                      stop("Input (tp1) is of wrong class.")
                    }
                    if (!is(group1, "character")) {
                      stop("Input (group1) is of wrong class.")
                    }
                    if ((ncol(data1) != length(tp1)) & 
                        (ncol(data1) != length(group1))) {
                      stop("Dimension do not correspond")
                    }
                    if (!is(data2, "matrix")) {
                      stop("Input (data2) is of wrong class.")
                    }
                    if (!is(tp2, "factor")) {
                      stop("Input (tp2) is of wrong class.")
                    }
                    if (!is(group2, "character")) {
                      stop("Input (group2) is of wrong class.")
                    }
                    if ((ncol(data2) != length(tp2)) & 
                        (ncol(data2) != length(group2))) {
                      stop("Dimension of the design matrices do not correspond")
                    }
                    plot_panel <- function(data,
                                           tp,
                                           diet,
                                           gene){
                      
                      d <- tibble(time=tp, 
                                  value=as.numeric(data[rownames(data) == gene_name, ]),
                                  diet=diet
                      ) %>% arrange(time)
                      
                      df <- data_summary(d,
                                         varname="value", 
                                         groupnames=c("diet", "time")
                      )
                      
                      p1 <- ggplot(df, aes(x=time, 
                                           y=mean,
                                           group=diet, 
                                           color=diet)) + 
                        geom_line() + 
                        geom_point() + 
                        ggtitle("hippocampus") + 
                        xlab("time points") + 
                        ylab("expression") +
                        theme(plot.title=element_text(hjust=0.5)) +
                        geom_errorbar(aes(ymin=mean - sd, 
                                          ymax=mean + sd), 
                                      width=.2, 
                                      position=position_dodge(0.05)) + 
                        scale_color_manual(values=c("lightskyblue", 
                                                      "indianred1")) +
                        annotate("rect", 
                                 xmin=3.5, 
                                 xmax=Inf, 
                                 ymin=-Inf,
                                 ymax=Inf, 
                                 alpha=.3
                        )
                    }
                    
                    p1 <- plot_panel(data=data1,
                                     tp=tp1,
                                     diet=group1,
                                     gene=gene_name
                    )
                    p2 <- plot_panel(data=data2,
                                     tp=tp2,
                                     diet=group2,
                                     gene=gene_name
                    )
                    
                    grid.arrange(p1,# + ylim(min(hip_data[rownames(hip_data)==gene,])-0.2, max(hip_data[rownames(hip_data)==gene,])+0.2), 
                                 p2,# + ylim(min(hpt_data[rownames(hpt_data)==gene,])-0.2, max(hpt_data[rownames(hpt_data)==gene,])+0.2),
                                 nrow=1, 
                                 top=gene_name
                    )
}
