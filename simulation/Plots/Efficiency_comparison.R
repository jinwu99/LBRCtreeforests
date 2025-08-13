rm(list=ls())

setwd("C:/~")

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(stringr)
library(emmeans)
library(patchwork)
library(broom)

method_names <- c(
  "LTRCctree",
  "LBRCctreeCC",
  "LBRCctreeCF",
  "LBRCctreeFF",
  "LBRCctreeFC"
  )
method_lookup <- data.frame(Method = method_names, stringsAsFactors = FALSE) %>%
  mutate(
    MethodLabel = factor(
      x = c(
        "LTRCctree",
        "LBRCctreeCC",
        "LBRCctreeCF",
        "LBRCctreeFF",
        "LBRCctreeFC"
      ),
      levels = c(
        "LTRCctree",
        "LBRCctreeCC",
        "LBRCctreeCF",
        "LBRCctreeFF",
        "LBRCctreeFC"
      )
    )
  )


model_plot <- function(models){
  sim <- 500
  dummy_L2 <- list()
  if(models=='tree'){
    distributions <- c("WI", "WD", "Lgn","Bat")
  }else{
    distributions <- c("WI", "WD")
  }

  Ns <- c(100,200,400)
  for(model in models) {
    for(dist in distributions) {
      for(N in Ns) {
        for(C in c(20,50)) {
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
            dummy_L2[[length(dummy_L2)+1]] <- data.frame(
              Model        = model,
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

  df_L2 <- bind_rows(dummy_L2)
  df_L2 <- df_L2[complete.cases(df_L2),]

  df_L2 <- df_L2 %>%
    left_join(method_lookup, by = "Method") %>%
    mutate(
      perm_method = factor(str_sub(MethodLabel, 10, 10)),
      pred_method = factor(str_sub(MethodLabel, 11, 11))
    )

  df_L2 <- df_L2 %>%
    group_by(Model, Distribution, N, C, MethodLabel) %>%
    mutate(rep = row_number()) %>%
    ungroup()

  baseline <- df_L2 %>%
    filter(MethodLabel == "LTRCctree") %>%
    select(Model, Distribution, N, C, rep, baseline_L2 = L2)

  df_long <- df_L2 %>%
    left_join(baseline,
              by = c("Model","Distribution","N","C","rep")) %>%
    mutate(diff = (baseline_L2-L2) )

  df_long <- df_long %>%
    mutate(across(c(Model, Distribution, N, C),
                  ~ factor(.x,
                           levels = if(is.numeric(.x)) sort(unique(.x)) else unique(.x)
                  )))

  df_long <- df_long[df_long$Method != "LTRCctree",]

  # fit the ANOVA (Distribution already filtered out for subsets;
  anova_fit <- aov(diff ~ Distribution + N + C + perm_method + pred_method,
                   data = df_long)
  print(summary(anova_fit))

  factors <- c(
    # "Distribution",
    "N","C","perm_method","pred_method")
  emm_list <- setNames(
    lapply(factors, function(fac) emmeans(anova_fit, as.formula(paste0("~", fac)))),
    factors
  )

  plot_df <- bind_rows(
    lapply(names(emm_list), function(fac) {
      as.data.frame(emm_list[[fac]]) %>%
        rename(Level = all_of(fac)) %>%
        mutate(Factor = fac)
    })
  ) %>%
    mutate(Factor = recode(Factor,
                           # Distribution = "Distribution",
                           N            = "Sample size",
                           C            = "Censoring rate",
                           perm_method  = "Variable selection",
                           pred_method  = "Survival prediction"
    )) %>%
    mutate(Factor = factor(Factor,
                           levels = c(
                                      # "Distribution",
                                      "Sample size",
                                      "Censoring rate",
                                      "Variable selection",
                                      "Survival prediction")
    ))

  ref_lines <- plot_df %>%
    group_by(Factor) %>%
    summarise(grand = mean(emmean, na.rm = TRUE), .groups = "drop")

  title <- paste0(toupper(substring(models, 1, 1)), substring(models, 2))
  plot_df$title <- title

  gg <- ggplot(plot_df, aes(x = Level, y = emmean, group = 1)) +
    geom_point(size = 2, color = "steelblue") +
    geom_line(color = "steelblue") +
    geom_hline(data = ref_lines,
               aes(yintercept = grand),
               linetype = "dashed", color = "grey40") +
    facet_grid(title ~ Factor,
               scales = "free_x",
               space = "free_x") +
    labs(
      x     = NULL,
      y     = "Mean"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
      strip.background = element_rect(color = "black", fill = "white", size = 0.5),
      strip.text       = element_text(size = 14),
      axis.text.x      = element_text(hjust = 1, size=13),
      axis.text.y      = element_text(size=12),
      panel.spacing.x  = unit(0.3, "lines"),
      axis.title.y = element_text(size = 13)
    )

  gg <- gg + geom_hline(yintercept = 0, color = "black", size = 0.5)

  anova_tbl <- broom::tidy(anova_fit) %>%
    filter(term %in% c("N","C","perm_method","pred_method")) %>%
    select(term, sumsq, p.value) %>%
    mutate(
      sumsq = round(sumsq, 4),
      p     = ifelse(p.value < 0.001, "< 0.001", paste0("= ", round(p.value, 3)))
    )

  annotation_df <- anova_tbl %>%
    mutate(
      Factor = recode(term,
                      N           = "Sample size",
                      C           = "Censoring rate",
                      perm_method = "Variable selection",
                      pred_method = "Survival prediction"
      ),
      Factor = factor(Factor, levels = levels(plot_df$Factor)),
      label = paste0("Sum Sq = ", sumsq, "\n", "Pr(>F) ", p)
    )

  gg <- gg +
    geom_text(
      data      = annotation_df,
      aes(x = -Inf, y = -Inf, label = label),
      hjust     = -0.1,
      vjust     = -0.5,
      size      = 3.5,
      inherit.aes = FALSE
    )


  return(gg)
}

ggtree <- model_plot("tree")
gglinear <- model_plot("linear")
ggnonlinear <- model_plot("nonlinear")

final_plot <- ggtree / gglinear / ggnonlinear +
  plot_layout(heights = c(1, 1, 1))

final_plot



