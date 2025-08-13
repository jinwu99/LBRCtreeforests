rm(list=ls())

setwd("C:/~")

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(simplecolors)

models <- c("tree","linear","nonlinear")
dists <- c("WI","WD")
# models <- c("tree")
# dists <- c("Lgn","Bat")
Ns <- c(100,200,400)
Cs <- c(20)

sim <- 500
dummy_L2 <- list()

storage <- data.frame(
  model = character(),
  dist = character(),
  N = numeric(),
  C = numeric(),
  mtry = numeric(),
  method = character(),
  mean = numeric(),
  sd = numeric(),
  se = numeric()
)

th <- 500

for(model in models){
  for(dist in dists) {
    # if(dist %in% c("Lgn","Bat") & model %in% c("linear","nonlinear")) next
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

        tmp <- data.frame(model=NA, dist=NA, N=NA, C=NA, mtry=NA,
                          method=NA, mean=NA, sd=NA, se=NA)
        tmp$model <- model
        tmp$dist <- dist
        tmp$N <- N
        tmp$C <- C

        mtry_candidates <- c(1,2,3,6,12,24,30,0)
        for(method in c("LTRCcforest","LBRCcforestCC","LBRCcforestFF")){
          tmp$method <- method
          for(m in mtry_candidates){
            if(m==0){
              tmp_mtry <- "mtryOPT"
              tmp$mtry <- "Tuned"
            }else{
              tmp_mtry <- paste0("mtry",m)
              tmp$mtry <- paste(m)
            }

            tmp_method <- abs(RES[["L2"]][[method]][[tmp_mtry]][1:th] - RES[["L2"]][[method]][["mtryMIN"]][1:th])
            tmp$mean <- mean(tmp_method)
            tmp$sd <- sd(tmp_method)
            tmp$se <- tmp$sd/sqrt(th)
            storage <- rbind(storage,tmp)
          }
        }


      }
    }
  }
}

method_names <- c(
  "LTRCcforest", "LBRCcforestCC", "LBRCcforestFF"
)
method_lookup <- data.frame(method = method_names, stringsAsFactors = FALSE) %>%
  mutate(
    Model = factor(
      x = c(
        "LTRC-CIF",
        "LBRC-CIF-C",
        "LBRC-CIF-F"
      ),
      levels = c(
        "LTRC-CIF",
        "LBRC-CIF-C",
        "LBRC-CIF-F"
      )
    )
  )

mtrys <- c(paste(mtry_candidates),"Tuned")

plot_df <- storage %>%
  left_join(method_lookup, by = "method") %>%
  mutate(
    mtry = factor(mtry, levels = mtrys),
    facet_label = paste0(dist, ", ", model, ", n = ", N)
  )

desired_order <- c(
  "WI, tree, n = 100",
  "WI, tree, n = 200",
  "WI, tree, n = 400",
  "WD, tree, n = 100",
  "WD, tree, n = 200",
  "WD, tree, n = 400",
  "WI, linear, n = 100",
  "WI, linear, n = 200",
  "WI, linear, n = 400",
  "WD, linear, n = 100",
  "WD, linear, n = 200",
  "WD, linear, n = 400",
  "WI, nonlinear, n = 100",
  "WI, nonlinear, n = 200",
  "WI, nonlinear, n = 400",
  "WD, nonlinear, n = 100",
  "WD, nonlinear, n = 200",
  "WD, nonlinear, n = 400"
)

# desired_order <- c(
#   "Lgn, tree, n = 100",
#   "Lgn, tree, n = 200",
#   "Lgn, tree, n = 400",
#   "Bat, tree, n = 100",
#   "Bat, tree, n = 200",
#   "Bat, tree, n = 400"
# )

plot_df <- plot_df %>%
  mutate(facet_label = factor(facet_label, levels = desired_order))

ggplot(plot_df, aes(x = mtry, y = mean, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  facet_wrap(~ facet_label,
             nrow = 3,
             ncol = 6,
             scales = "free_y") +
  # facet_wrap(~ facet_label,
  #            nrow = 2,
  #            ncol = 3,
  #            scales = "free_y") +
  scale_fill_manual(
    values = c(
      "LTRC-CIF"   = "olivedrab3",
      "LBRC-CIF-C"  = "dodgerblue",
      "LBRC-CIF-F"   = "tomato"
    )
  ) +
  # labels & theme
  labs(
    x = NULL,
    y = "Mean integrated L2 distance difference",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    legend.text     = element_text(size = 12),
    axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1, size=13),
    strip.text       = element_text(size = 14),
    axis.title.y     = element_text(size = 15, vjust=1),
    legend.position = "bottom",
    strip.background = element_rect(color = "black", fill = "white", size = 0.5),
  )















