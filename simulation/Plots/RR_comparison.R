rm(list=ls())

setwd("C:/~")

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)

model <- "tree"
sim <- 500
dummy_data <- list()
distributions <- c("WI", "WD", "Lgn", "Bat")
for(dist in distributions) {
  for(N in c(100,200,400)) {
    for(C in c(20,50)) {
      P_vec <- c(30)
      scores <- list()
      scores$LTRCctree <- c()
      scores$LBRCctreeH <- c()
      scores$LBRCctreeF <- c()
      for(m in 1:length(P_vec)){
        P <- P_vec[m]
        base_location <- "./results"
        final_location <- paste0(base_location,"/",model,"/",dist,"/N",N)
        setwd(final_location)
        string = paste0("LBRC_DIST_",dist,
                        "_MODEL_",model,
                        "_P",P,
                        "_N",N,
                        "_C",C,
                        "_sim",sim)
        load(file=string)
        scores$LTRCctree[m] <- mean(RES$RR$LTRCctree,na.rm=T)
        scores$LBRCctreeC[m] <- mean(RES$RR$LBRCctreeC,na.rm=T)
        scores$LBRCctreeF[m] <- mean(RES$RR$LBRCctreeF,na.rm=T)
      }
      dummy_data[[length(dummy_data) + 1]] <- list(
        Distribution = dist,
        N = N,
        C = C,
        LTRCctree   = scores$LTRCctree,
        LBRCctreeC   = scores$LBRCctreeC,
        LBRCctreeF  = scores$LBRCctreeF
      )
    }
  }
}

df <- do.call(rbind, lapply(dummy_data, function(x) {
  data.frame(
    Distribution = x$Distribution,
    N = x$N,
    C = x$C,
    P = c(30),
    LTRCctree = x$LTRCctree,
    LBRCctreeC = x$LBRCctreeC,
    LBRCctreeF = x$LBRCctreeF
  )
}))

colnames(df)

colnames(df)[5] <- "LTRC-CIT"
colnames(df)[6] <- "LBRC-CIT(C,·)"
colnames(df)[7] <- "LBRC-CIT(F,·)"

df_long <- df %>%
  pivot_longer(cols = c(
                        "LTRC-CIT",
                        "LBRC-CIT(C,·)",
                        "LBRC-CIT(F,·)"),
               names_to = "Model",
               values_to = "RR")


df_long$Model <- factor(
  df_long$Model,
  levels = c(
    "LTRC-CIT",
    "LBRC-CIT(C,·)",
    "LBRC-CIT(F,·)"
  )
)



big_RR_plt2 <- function(dist){
  df_exp <- df_long %>% filter(Distribution == dist) %>% mutate(
    facet_label = paste0(Distribution,", c = ", round(C/100,1))
  )

  plot_by_nc <- function(data) {
    p <- ggplot(data, aes(x = factor(N), y = RR, group = Model, color = Model)) +
      geom_line(
        data     = filter(data, Model == "LTRC-CIT"),
        aes(x = as.numeric(factor(N)), y = RR),
        alpha    = 0.7
      ) +
      geom_line(
        data     = filter(data, Model == "LBRC-CIT(C,·)"),
        aes(x = as.numeric(factor(N)) + -0.05, y = RR),
        alpha    = 0.7
      ) +
      geom_line(
        data     = filter(data, Model == "LBRC-CIT(F,·)"),
        aes(x = as.numeric(factor(N)) +  0.05, y = RR),
        alpha    = 0.7
      ) +
      geom_point(
        data     = filter(data, Model == "LTRC-CIT"),
        aes(x = as.numeric(factor(N)), y = RR),
        alpha    = 0.7
      ) +
      geom_point(
        data     = filter(data, Model == "LBRC-CIT(C,·)"),
        aes(x = as.numeric(factor(N)) + -0.05, y = RR),
        alpha    = 0.7
      ) +
      geom_point(
        data     = filter(data, Model == "LBRC-CIT(F,·)"),
        aes(x = as.numeric(factor(N)) +  0.05, y = RR),
        alpha    = 0.7
      ) +
      coord_cartesian(ylim = c(0:1)) +
      facet_wrap(
        ~ facet_label,
        ncol = 1
      ) +
      scale_x_discrete(
        limits = as.character(sort(unique(data$N))),
        labels =            sort(unique(data$N))
      ) +
      scale_y_continuous(
        breaks = scales::pretty_breaks(),
        labels = scales::number_format(accuracy = 0.1)
      ) +
      scale_color_manual(
        values = c(
          "LTRC-CIT"   = "olivedrab3",
          "LBRC-CIT(C,·)"  = "dodgerblue",
          "LBRC-CIT(F,·)"   = "tomato"
        )
      ) +
      labs(
        x = "Sample Size",
        y = NULL
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 0.5),
        strip.text       = element_text(size = 14),
        axis.text.x      = element_text(hjust = 0.5, size=11),
        axis.text.y      = element_text(size=12),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 13),
        legend.text     = element_text(size = 12),

      )
    return(p)
  }

  plot_by_nc(df_exp)
}


p_list <- list()
for(dist in c("WD","WI","Lgn","Bat")){
  p_list[[dist]] <- big_RR_plt2(dist)
}

all_dists <- ggarrange(
  p_list$WI, p_list$WD, p_list$Lgn, p_list$Bat,
  ncol = 4
)

legend_plot <- ggplot(df_long, aes(x = factor(N), y = RR, group = Model, color = Model)) +
  geom_line() +
  geom_point() +
  scale_color_manual(
    values = c(
      "LTRC-CIT"   = "olivedrab3",
      "LBRC-CIT(C,·)"  = "dodgerblue",
      "LBRC-CIT(F,·)"   = "tomato"
    )
  ) +
  theme_bw() +
  theme(
    legend.position      = "bottom",
    legend.text          = element_text(size = 13),
    legend.margin        = margin(t = 0, b = 0, unit = "pt"),
    legend.key.height    = unit(3, "pt"),
    legend.spacing.y     = unit(0, "pt")
  ) +
  guides(
    color = guide_legend(
      keyheight   = unit(6, "pt"),
      default.unit = "pt"
    )
  )
common_legend <- get_legend(legend_plot)

all_dists_with_legend <- ggarrange(
  all_dists,
  common_legend,
  ncol = 1,
  heights = c(1, 0.05)
)

final_plot <- annotate_figure(
  all_dists_with_legend,
  left = text_grob("Recovery Rate", rot = 90, size = 15)
)

final_plot


