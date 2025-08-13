rm(list=ls())

setwd("C:/~")

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(simplecolors)

method_names <- c(
  "KM_LT",  "MCLE", "MFLE",
  "LTRCctree", "LBRCctreeCC", "LBRCctreeFF",
  "LTRCcforest", "LBRCcforestCC", "LBRCcforestFF",
  "LTRCcox",
  "LBRCcox"
)
method_lookup <- data.frame(Method = method_names, stringsAsFactors = FALSE) %>%
  mutate(
    MethodLabel = factor(
      x = c(
        "L1",
        "C1",
        "F1",
        "LT",
        "CT",
        "FT",
        "LF",
        "CF",
        "FF",
        "LC",
        "BC"
      ),
      levels = c(
        "L1",
        "C1",
        "F1",
        "LT",
        "CT",
        "FT",
        "LF",
        "CF",
        "FF",
        "LC",
        "BC"
      )
    )
  )


L2_plot <- function(models,dists,Ns,Cs,num_rc){
  sim <- 500
  dummy_L2 <- list()
  for(model in models){
    for(dist in dists) {
      for(N in Ns){
        for(C in Cs) {
          base_location <- "./results"
          final_location <- paste0(base_location,"/",model,"/",dist,"/N",N)
          setwd(final_location)
          fname <- paste0("LBRC_DIST_",dist,
                          "_MODEL_",model,
                          "_P30",
                          "_N",N,
                          "_C",C,
                          "_sim",sim)
          load(fname)

          for(m in method_names) {
            if(m %in% c("LTRCcforest","LBRCcforestCC","LBRCcforestFF")){
              dummy_L2[[length(dummy_L2)+1]] <- data.frame(
                Model        = paste0(toupper(substring(model, 1, 1)), substring(model, 2)),
                Distribution = dist,
                N            = N,
                C            = C,
                Method       = m,
                L2           = RES$L2[[m]]$mtryOPT
              )
            }else{
              dummy_L2[[length(dummy_L2)+1]] <- data.frame(
                Model        = paste0(toupper(substring(model, 1, 1)), substring(model, 2)),
                Distribution = dist,
                N            = N,
                C            = C,
                Method       = m,
                L2           = RES$L2[[m]]
              )
            }
          }
        }
      }
    }
  }

  df_L2 <- bind_rows(dummy_L2)
  df_L2 <- df_L2[complete.cases(df_L2),]
  df_L2$log_L2 <- log(df_L2$L2)

  df_L2 <- df_L2 %>%
    left_join(method_lookup, by = "Method") %>%
    mutate(
      facet_label = paste0(Distribution,", ",Model, ", (", N,", ", round(C/100,1),")")
    )

  df_L2$facet_label <- factor(df_L2$facet_label,
                              levels = unique(df_L2$facet_label))

  p <- ggplot(df_L2, aes(x = MethodLabel, y = log_L2, fill = MethodLabel)) +
    geom_boxplot(outlier.size = 0.8) +
    scale_fill_manual(
      values = sc(
        "L1" = "grey1",  "C1" = "grey2",  "F1" = "grey3",
        "LT" = "green1", "CT" = "green2", "FT" = "green3",
        "LF" = "blue1",  "CF" = "blue2",  "FF" = "blue3",
        "LC" = "yellow1","BC" = "yellow2"
      )
    ) +
    facet_wrap(~ facet_label, nrow = num_rc[1], ncol=num_rc[2], scales = "free_y") +
    scale_y_continuous(
      breaks = scales::pretty_breaks(),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    labs(
      x = "Method",
      y = "Log of integrated L2 distance"
    ) +
    theme_bw() +
    theme(
      strip.text       = element_text(size = 14),
      axis.text.x      = element_text(size = 13, angle=90, hjust=0, vjust=0.5),
      axis.title.y     = element_text(size = 15, vjust=1),
      axis.title.x     = element_text(size = 15, vjust=-1),
      legend.position  = "none",
      strip.background = element_rect(color = "black", fill = "white", size = 0.5),
    )
  return(p)
}

L2_plot(c("tree","linear","nonlinear"),c("WI"),c(100,200,400),c(20,50),c(4,6))
# L2_plot(c("tree","linear","nonlinear"),c("WD"),c(100,200,400),c(20,50),c(6,6))
# L2_plot(c("tree"),c("Lgn","Bat"),c(100,200,400),c(20,50),c(6,6))
