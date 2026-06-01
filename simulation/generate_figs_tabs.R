library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(simplecolors)
library(stringr)
library(emmeans)
library(patchwork)
library(broom)
library(scales)

##### figure 1: recovery rate comparison across tree models #####
RR_plot <- function(results_dir){
  model <- "tree"
  distributions <- c("WI", "WD", "Lgn", "Bat")

  dummy_data <- list()
  for(dist in distributions) {
    for(N in c(100,200,400)) {
      for(C in c(20,50)) {
        scores <- list()
        scores$LTRCctree <- c()
        scores$LBRCctreeH <- c()
        scores$LBRCctreeF <- c()

        result_location <- paste0(results_dir,"/",model,"/",dist,"/N",N)
        setwd(result_location)
        fname <- paste0("LBRC_DIST_",dist,
                        "_MODEL_",model,
                        "_P30",
                        "_N",N,
                        "_C",C)
        load(fname)

        scores$LTRCctree <- mean(RES$RR$LTRCctree,na.rm=T)
        scores$LBRCctreeC <- mean(RES$RR$LBRCctreeC,na.rm=T)
        scores$LBRCctreeF <- mean(RES$RR$LBRCctreeF,na.rm=T)
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
      P = 30,
      LTRCctree = x$LTRCctree,
      LBRCctreeC = x$LBRCctreeC,
      LBRCctreeF = x$LBRCctreeF
    )
  }))

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
          alpha    = 1,
          linewidth = 1.2
        ) +
        geom_line(
          data     = filter(data, Model == "LBRC-CIT(C,·)"),
          aes(x = as.numeric(factor(N)) + -0.05, y = RR),
          alpha    = 0.7,
          linewidth = 1.2
        ) +
        geom_line(
          data     = filter(data, Model == "LBRC-CIT(F,·)"),
          aes(x = as.numeric(factor(N)) +  0.05, y = RR),
          alpha    = 0.7,
          linewidth = 1.2
        ) +
        geom_point(
          data     = filter(data, Model == "LTRC-CIT"),
          aes(x = as.numeric(factor(N)), y = RR),
          alpha    = 1,
          size = 3
        ) +
        geom_point(
          data     = filter(data, Model == "LBRC-CIT(C,·)"),
          aes(x = as.numeric(factor(N)) + -0.05, y = RR),
          alpha    = 0.7,
          size = 3
        ) +
        geom_point(
          data     = filter(data, Model == "LBRC-CIT(F,·)"),
          aes(x = as.numeric(factor(N)) +  0.05, y = RR),
          alpha    = 0.7,
          size = 3
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
          panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
          strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
          strip.text       = element_text(size = 17),
          axis.text.x      = element_text(size=14),
          axis.text.y      = element_text(size=14),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 15),
          legend.text     = element_text(size = 14)
        )
      return(p)
    }

    return(plot_by_nc(df_exp))
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
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c(
        "LTRC-CIT"   = "olivedrab3",
        "LBRC-CIT(C,·)"  = "dodgerblue",
        "LBRC-CIT(F,·)"   = "tomato"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.margin = margin(t = 0, b = 0, unit = "pt"),
      legend.key.width = unit(25, "pt")
    ) +

    guides(
      color = guide_legend(
        keyheight   = unit(6, "pt"),
        default.unit = "pt"
      )
    )

  all_dists_with_legend <- ggarrange(
    all_dists,
    get_legend(legend_plot),
    ncol = 1,
    heights = c(1, 0.05)
  )

  p <- annotate_figure(
    all_dists_with_legend,
    left = text_grob("Recovery Rate", rot = 90, size = 18)
  )

  setwd(results_dir)
  ggsave(filename = "figure_1_compare_recovery_rate.pdf",
         plot = p, width = 12, height = 6, units = "in")

  cat("figure_1_compare_recovery_rate.pdf generated\n")
}


##### figure 2: efficiency gain of LBRCtrees over LTRCtree #####
ANOVA_plot <- function(results_dir){
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

  anova_model_plot <- function(models){
    sim <- 500
    dummy_L2 <- list()
    if(models=='tree'){
      distributions <- c("WI", "WD", "Lgn","Bat")
    }else{
      distributions <- c("WI", "WD")
    }

    for(model in models) {
      for(dist in distributions) {
        for(N in c(100,200,400)) {
          for(C in c(20,50)) {
            result_location <- paste0(results_dir,"/ANOVA_study/",model,"/",dist,"/N",N)
            setwd(result_location)
            fname <- paste0("LBRC_DIST_",dist,
                            "_MODEL_",model,
                            "_P30",
                            "_N",N,
                            "_C",C)
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
      mutate(diff = ((baseline_L2-L2)/baseline_L2) )
      # mutate(diff = (baseline_L2-L2) )

    df_long <- df_long %>%
      mutate(across(c(Model, Distribution, N, C),
                    ~ factor(.x,
                             levels = if(is.numeric(.x)) sort(unique(.x)) else unique(.x)
                    )))

    df_long <- df_long[df_long$Method != "LTRCctree",]

    # fit the ANOVA (Distribution already filtered out for subsets)
    anova_fit <- aov(diff ~ Distribution + N + C + perm_method + pred_method,
                     data = df_long)

    factors <- c("N","C","perm_method","pred_method")
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
                             N            = "Sample size",
                             C            = "Censoring rate",
                             perm_method  = "Variable selection",
                             pred_method  = "Survival prediction"
      )) %>%
      mutate(Factor = factor(Factor,
                             levels = c(
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
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                    width = 0.2, color = "steelblue", linewidth = 0.8) +
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
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        strip.text       = element_text(size = 15),
        axis.text.x      = element_text(hjust = 1, size=14),
        axis.text.y      = element_text(size = 13),
        panel.spacing.x  = unit(0.3, "lines"),
        axis.title.y = element_text(size = 15)
      )

    gg <- gg + geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

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
        size      = 4.25,
        inherit.aes = FALSE
      )

    return(gg)
  }

  ggtree <- anova_model_plot("tree")
  gglinear <- anova_model_plot("linear")
  ggnonlinear <- anova_model_plot("nonlinear")
  gginteraction <- anova_model_plot("interaction")

  p <- ggtree / gglinear / ggnonlinear / gginteraction +
    plot_layout(heights = c(1, 1, 1, 1))

  setwd(results_dir)
  ggsave(filename = paste0(mode,".pdf"),
         plot = p, width = 12, height = 11, units = "in")

  cat(paste0(mode,".pdf"), "generated\n")
}


##### figure 3 & S2 ~ S6: validation of OOB tuning for LBRCforests #####
OOBtune_plot <- function(mode, tune.metric = "brier", results_dir){
  if(mode %in% c("figure_3_WI_LBRCforests_OOB_tuning_brier",
                 "figure_S2_WI_LBRCforests_OOB_tuning_cindex")){
    dists <- c("WI")
    models <- c("tree","linear","nonlinear","interaction")
    width <- 12; height <- 8
    nrow <- 3; ncol <- 4
    facet_dir <- "v"
  }else if(mode %in% c("figure_S3_WD_LBRCforests_OOB_tuning_brier",
                       "figure_S4_WD_LBRCforests_OOB_tuning_cindex")){
    dists <- c("WD")
    models <- c("tree","linear","nonlinear","interaction")
    width <- 12; height <- 8
    nrow <- 3; ncol <- 4
    facet_dir <- "v"
  }else if(mode %in% c("figure_S5_LgnBat_LBRCforests_OOB_tuning_brier",
                       "figure_S6_LgnBat_LBRCforests_OOB_tuning_cindex")){
    dists <- c("Lgn","Bat")
    models <- c("tree")
    width <- 10; height <- 6
    nrow <- 2; ncol <- 3
    facet_dir <- "h"
  }

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

  for(model in models){
    for(dist in dists) {
      for(N in c(100,200,400)){
        result_location <- paste0(results_dir,"/",model,"/",dist,"/N",N)
        setwd(result_location)

        fname <- paste0("LBRC_DIST_",dist,
                        "_MODEL_",model,
                        "_P30",
                        "_N",N,
                        "_C",20)
        load(fname)

        tmp <- data.frame(model=NA, dist=NA, N=NA, mtry=NA,
                          method=NA, mean=NA, sd=NA, se=NA)
        tmp$model <- paste0(toupper(substring(model, 1, 1)), substring(model, 2))
        tmp$dist <- dist
        tmp$N <- N

        mtry_candidates <- c(1,2,3,6,12,24,30,0)
        for(method in c("LTRCcforest","LBRCcforestCC","LBRCcforestFF")){
          tmp$method <- method
          for(m in mtry_candidates){
            if(m==0){
              tmp_mtry <- ifelse(tune.metric=="brier","mtryOPT","mtryOPT2")
              mtry_label <- "Tuned"
              tmp$mtry <- "Tuned"
            }else{
              tmp_mtry <- paste0("mtry",m)
              tmp$mtry <- paste(m)
            }


            tmp_method <- abs(RES[["L2"]][[method]][[tmp_mtry]] - RES[["L2"]][[method]][["mtryMIN"]] )
            tmp$mean   <- mean(tmp_method, na.rm=T)
            tmp$sd     <- sd(tmp_method, na.rm=T)
            tmp$se     <- tmp$sd / sqrt(sum(!is.na(tmp_method)))
            storage    <- rbind(storage, tmp)
          }
        }
      }
    }
  }


  method_names <- c("LTRCcforest", "LBRCcforestCC", "LBRCcforestFF")
  method_lookup <- data.frame(method = method_names, stringsAsFactors = FALSE) %>%
    mutate(
      Model = factor(
        x      = c("LTRC-CIF", "LBRC-CIF-C", "LBRC-CIF-F"),
        levels = c("LTRC-CIF", "LBRC-CIF-C", "LBRC-CIF-F")
      )
    )

  mtrys <- c(paste(mtry_candidates), "Tuned")

  plot_df <- storage %>%
    left_join(method_lookup, by = "method") %>%
    mutate(
      mtry = factor(mtry, levels = mtrys),
      facet_label = paste0(dist, ", ", model, ", n = ", N)
    )

  if(mode %in% c("figure_S5_LgnBat_LBRCforests_OOB_tuning",
                 "figure_S6_LgnBat_LBRCforests_OOB_tuning_cindex")){
    desired_order <- c(
      "Lgn, Tree, n = 100",
      "Lgn, Tree, n = 200",
      "Lgn, Tree, n = 400",
      "Bat, Tree, n = 100",
      "Bat, Tree, n = 200",
      "Bat, Tree, n = 400"
    )
    plot_df <- plot_df %>%
      mutate(facet_label = factor(facet_label, levels = desired_order))
  }else{
    plot_df <- plot_df %>%
      mutate(facet_label = factor(facet_label, levels = unique(plot_df$facet_label)))
  }


  p <- ggplot(plot_df, aes(x = mtry, y = mean, fill = Model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    facet_wrap(~ facet_label, nrow = nrow, ncol = ncol, scales = "free_y", dir = facet_dir) +
    scale_fill_manual(
      values = c(
        "LTRC-CIF"   = "olivedrab3",
        "LBRC-CIF-C" = "dodgerblue",
        "LBRC-CIF-F" = "tomato"
      )
    ) +
    labs(
      x = "mtry value",
      y = "Mean integrated L2 distance difference",
      fill = "Model"
    ) +
    theme_bw() +
    theme(
      legend.title     = element_text(size = 15),
      legend.text      = element_text(size = 15),
      axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
      axis.text.y      = element_text(size = 12),
      strip.text       = element_text(size = 15),
      axis.title.x     = element_text(size = 15, vjust = -1),
      axis.title.y     = element_text(size = 18),
      legend.position  = "bottom",
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    )

  setwd(results_dir)
  ggsave(filename = paste0(mode,".pdf"),
         plot = p, width = width, height = height, units = "in")

  cat(paste0(mode,".pdf"), "generated\n")
}


##### figure 4 & S7 ~ S10: prediction comparison across models #####
L2_plot <- function(mode, results_dir){
  if(mode == "figure_4_WI_20_compare_prediction_accuracy"){ # WI, C = 20%
    dists <- "WI"
    models <- c("tree","linear","nonlinear","interaction")
    cens_rates <- 20
  }else if(mode == "figure_S7_WI_50_compare_prediction_accuracy"){ # WI, C = 50%
    dists <- "WI"
    models <- c("tree","linear","nonlinear","interaction")
    cens_rates <- 50
  }else if(mode == "figure_S8_WD_20_compare_prediction_accuracy"){ # WD, C = 20%
    dists <- "WD"
    models <- c("tree","linear","nonlinear","interaction")
    cens_rates <- 20
  }else if(mode == "figure_S9_WD_50_compare_prediction_accuracy"){ # WD, C = 50%
    dists <- "WD"
    models <- c("tree","linear","nonlinear","interaction")
    cens_rates <- 50
  }else if(mode == "figure_S10_LgnBat_2050_compare_prediction_accuracy"){ # Bat, Lgn, C = 20%, 50%
    dists <- c("Lgn", "Bat")
    models <- "tree"
    cens_rates <- c(20, 50)
  }

  method_names <- c(
    "KM_LT",  "MCLE", "MFLE",
    "LTRCctree", "LBRCctreeCC", "LBRCctreeFF",
    "LTRCcforest", "LBRCcforestCC", "LBRCcforestFF",
    "LTRCcox", "LBRCcox"
  )
  method_lookup <- data.frame(Method = method_names, stringsAsFactors = FALSE) %>%
    mutate(
      MethodLabel = factor(
        x = c(
          "L1","C1","F1",
          "LT","CT","FT",
          "LF","CF","FF",
          "LC","BC"
        ),
        levels = c(
          "L1","C1","F1",
          "LT","CT","FT",
          "LF","CF","FF",
          "LC","BC"
        )
      )
    )

  dummy_L2 <- list()
  for(model in models){
    for(dist in dists){
      for(N in c(100,200,400)){
        for(C in cens_rates) {
          result_location <- paste0(results_dir,"/",model,"/",dist,"/N",N)
          setwd(result_location)
          fname <- paste0("LBRC_DIST_",dist,
                          "_MODEL_",model,
                          "_P30",
                          "_N",N,
                          "_C",C)
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

  df_L2$facet_label <- factor(df_L2$facet_label, levels = unique(df_L2$facet_label))
  if(mode == "figure_S10_LgnBat_2050_compare_prediction_accuracy"){
    desired_order <- c(
      "Lgn, Tree, (100, 0.2)", "Lgn, Tree, (200, 0.2)", "Lgn, Tree, (400, 0.2)",
      "Lgn, Tree, (100, 0.5)", "Lgn, Tree, (200, 0.5)", "Lgn, Tree, (400, 0.5)",
      "Bat, Tree, (100, 0.2)", "Bat, Tree, (200, 0.2)", "Bat, Tree, (400, 0.2)",
      "Bat, Tree, (100, 0.5)", "Bat, Tree, (200, 0.5)", "Bat, Tree, (400, 0.5)"
    )
    df_L2$facet_label <- factor(df_L2$facet_label, levels = desired_order)
  }else{
    df_L2$facet_label <- factor(df_L2$facet_label, levels = unique(df_L2$facet_label))
  }

  whisker_data <- df_L2 %>%
    group_by(facet_label, Method) %>%
    summarize(
      Q1 = quantile(log_L2, 0.25, na.rm = TRUE),
      Q3 = quantile(log_L2, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      upper_whisker = Q3 + 1.5 * IQR,
      .groups = "drop"
    ) %>%
    group_by(facet_label) %>%
    summarize(
      max_whisker = max(upper_whisker, na.rm = TRUE),
      .groups = "drop"
    )

  df_L2 <- df_L2 %>%
    left_join(whisker_data, by = "facet_label") %>%
    mutate(
      upper_limit = max_whisker + abs(max_whisker) * 0.05,
      log_L2_plot = ifelse(log_L2 > upper_limit, upper_limit, log_L2)
    )

  p <- ggplot(df_L2, aes(x = MethodLabel, y = log_L2_plot, fill = MethodLabel)) +
    geom_boxplot(
      outlier.size = 0.8,
      outlier.alpha = 0.4,
      outlier.color = "grey30"
    ) +
    scale_fill_manual(
      values = sc(
        "L1" = "grey1",  "C1" = "grey2",  "F1" = "grey3",
        "LT" = "green1", "CT" = "green2", "FT" = "green3",
        "LF" = "blue1",  "CF" = "blue2",  "FF" = "blue3",
        "LC" = "yellow1","BC" = "yellow2"
      )
    ) +
    facet_wrap(~ facet_label, nrow = 3, ncol = 4, scales = "free_y", dir = "v") +
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
      legend.position  = "none",
      axis.text.x      = element_text(size = 13, angle = 90, hjust = 0, vjust = 0.5),
      axis.text.y      = element_text(size = 15),
      axis.title.y     = element_text(size = 18),
      axis.title.x     = element_text(size = 17, vjust = -1),
      strip.text       = element_text(size = 14),
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
      plot.caption     = element_text(size = 10, color = "grey50")
    )

  setwd(results_dir)
  ggsave(filename = paste0(mode,".pdf"),
         plot = p, width = 12, height = 8, units = "in")

  cat(paste0(mode,".pdf"), "generated\n")
}


##### figure S1 & S11: test of unbiasedness of variable selection #####
Unbiasedness_table <- function(mode, results_dir){
  if(mode == "figure_S1_test_unbiasedness_LBRCtrees"){
    dist_codes  <- c("WI","WD","Lgn")
    dist_labels <- c(WI  = "Weibull-I",
                     WD  = "Weibull-D",
                     Lgn = "Lognormal")

    rows <- list()
    result_location <- paste0(results_dir,"/test_unbiasedness")
    setwd(result_location)

    for (C in c(20,50)) {
      for (dist in dist_codes) {
        fname <- paste0("LBRC_UNBIASED_TEST_DIST_",dist,
                        "_C",C)
        load(fname)

        # Ensure all X1–X6 appear in the tables (even if count = 0)
        ctable <- table(RES$RV$LBRCctreeC)
        ftable <- table(RES$RV$LBRCctreeF)
        cp <- chisq.test(ctable)$p.value
        fp <- chisq.test(ftable)$p.value

        rows[[length(rows) + 1]] <- data.frame(
          Censoring    = paste0(C, "%"),
          Distribution = dist_labels[dist],
          C_X1 = as.numeric(ctable[1]),
          C_X2 = as.numeric(ctable[2]),
          C_X3 = as.numeric(ctable[3]),
          C_X4 = as.numeric(ctable[4]),
          C_X5 = as.numeric(ctable[5]),
          C_X6 = as.numeric(ctable[6]),
          C_pval  = round(cp, 4),
          F_X1 = as.numeric(ftable[1]),
          F_X2 = as.numeric(ftable[2]),
          F_X3 = as.numeric(ftable[3]),
          F_X4 = as.numeric(ftable[4]),
          F_X5 = as.numeric(ftable[5]),
          F_X6 = as.numeric(ftable[6]),
          F_pval  = round(fp, 4),
          stringsAsFactors = FALSE
        )
      }
    }

    unbiased_tab <- do.call(rbind, rows)

    setwd(results_dir)
    csv_name <- "test_unbiasedness_LBRCtrees.csv"
    write.csv(unbiased_tab, file=csv_name, row.names=FALSE)

    cat("test_unbiasedness_LBRCtrees.csv generated\n")
  }else if(mode == "figure_S11_sensitivity_analysis_unbiasedness"){
    rows <- list()
    tau <- qweibull(0.9999,2,3) + 1
    # rho_grid = c(tau/10, tau/5, tau/2, tau/1)
    rho_grid = c(0, tau/10, tau/5, tau/2)

    for(scenario in c("unbias_texpt","unbias_covd")){
      for(rho in rho_grid){
        if(scenario == "unbias_texpt"){
          trunc_rate <- round(rho/tau,2)
          if(trunc_rate > 0){
            result_location <- paste0(results_dir,"/sensitivity_analysis/test_unbiasedness")
            setwd(result_location)
            string = sprintf("LBRC_UNBIASED_TEST_DIST_%s_C%1.0f_T%1.2f",
                             "WI",20,trunc_rate)
          }else{
            result_location <- paste0(results_dir,"/test_unbiasedness")
            setwd(result_location)
            string = sprintf("LBRC_UNBIASED_TEST_DIST_%s_C%1.0f",
                             "WI",20)
          }
        }else{ # scenario == "unbias_covd"
          trunc_rate <- round(rho/tau,2)
          if(trunc_rate > 0){
            result_location <- paste0(results_dir,"/sensitivity_analysis/test_unbiasedness")
            setwd(result_location)
            string = sprintf("LBRC_UNBIASED_TEST_DIST_%s_C%1.0f_S%1.2f",
                             "WI",20,trunc_rate)
          }else{
            result_location <- paste0(results_dir,"/test_unbiasedness")
            setwd(result_location)
            string = sprintf("LBRC_UNBIASED_TEST_DIST_%s_C%1.0f",
                             "WI",20)
          }
        }

        load(string)

        ltable <- table(RES$RV$LTRCctree)
        ctable <- table(RES$RV$LBRCctreeC)
        ftable <- table(RES$RV$LBRCctreeF)
        lp <- round(chisq.test(ltable)$p.value, 3)
        cp <- round(chisq.test(ctable)$p.value, 3)
        fp <- round(chisq.test(ftable)$p.value, 3)

        iteration_df <- data.frame(
          Scenario = rep(scenario, 3),
          mu       = rep(trunc_rate, 3),
          Method   = c("LTRCctree", "LBRCctreeC", "LBRCctreeF"), # Identifies the source
          X1       = c(as.numeric(ltable[1]), as.numeric(ctable[1]), as.numeric(ftable[1])),
          X2       = c(as.numeric(ltable[2]), as.numeric(ctable[2]), as.numeric(ftable[2])),
          X3       = c(as.numeric(ltable[3]), as.numeric(ctable[3]), as.numeric(ftable[3])),
          X4       = c(as.numeric(ltable[4]), as.numeric(ctable[4]), as.numeric(ftable[4])),
          X5       = c(as.numeric(ltable[5]), as.numeric(ctable[5]), as.numeric(ftable[5])),
          X6       = c(as.numeric(ltable[6]), as.numeric(ctable[6]), as.numeric(ftable[6])),
          pval     = c(lp, cp, fp),
          stringsAsFactors = FALSE
        )

        rows[[length(rows) + 1]] <- iteration_df
      }
    }

    unbiased_tab <- do.call(rbind, rows)

    setwd(results_dir)
    csv_name <- "sensitivity_unbiasedness.csv"
    write.csv(unbiased_tab, file=csv_name, row.names=FALSE)

    cat("sensitivity_unbiasedness.csv generated\n")
  }
}
Unbiasedness_figure <- function(mode, results_dir){
  fmt_p <- function(p){
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))
  }

  if(mode == "figure_S1_test_unbiasedness_LBRCtrees"){
    # results_dir <- paste0(results_dir,"/test_unbiasedness")
    dat <- read.csv(file.path(results_dir, "test_unbiasedness_LBRCtrees.csv"),
                    stringsAsFactors = FALSE)

    dat_long <- dat %>%
      mutate(
        Setting = paste0(
          case_when(
            Distribution == "Weibull-I" ~ "WI",
            Distribution == "Weibull-D" ~ "WD",
            Distribution == "Lognormal" ~ "Lgn",
            TRUE ~ Distribution
          ),
          ", ",
          gsub("%", "", Censoring), "%"
        )
      ) %>%
      select(Setting, starts_with("C_"), starts_with("F_")) %>%
      pivot_longer(
        cols = matches("^[CF]_X[1-6]$"),
        names_to = c("Method", "Variable"),
        names_pattern = "([CF])_(X[1-6])",
        values_to = "Frequency"
      ) %>%
      left_join(
        dat %>%
          mutate(
            Setting = paste0(
              case_when(
                Distribution == "Weibull-I" ~ "WI",
                Distribution == "Weibull-D" ~ "WD",
                Distribution == "Lognormal" ~ "Lgn",
                TRUE ~ Distribution
              ),
              ", ",
              gsub("%", "", Censoring), "%"
            )
          ) %>%
          select(Setting, C_pval, F_pval) %>%
          pivot_longer(
            cols = c(C_pval, F_pval),
            names_to = "Method",
            values_to = "pval"
          ) %>%
          mutate(Method = ifelse(Method == "C_pval", "C", "F")),
        by = c("Setting", "Method")
      ) %>%
      group_by(Setting, Method) %>%
      mutate(
        Total = sum(Frequency),
        Percent = Frequency / Total,
        Label = percent(Percent, accuracy = 0.1)
      ) %>%
      ungroup() %>%
      mutate(
        Setting = factor(
          Setting,
          levels = c("WI, 20%", "WD, 20%", "Lgn, 20%",
                     "WI, 50%", "WD, 50%", "Lgn, 50%")
        ),
        Method = factor(Method, levels = c("C", "F")),
        Variable = factor(Variable, levels = paste0("X", 1:6))
      )

    pdat <- dat_long %>%
      distinct(Setting, Method, Total, pval) %>%
      mutate(
        y = Total * 1.06,
        p_label = fmt_p(pval)
      )

    p <- ggplot(dat_long, aes(x = Method, y = Frequency, fill = Variable)) +
      geom_col(width = 0.75, color = "white", linewidth = 0.5) +
      geom_text(
        aes(label = Label),
        position = position_stack(vjust = 0.5),
        size = 5
      ) +
      geom_text(
        data = pdat,
        aes(x = Method, y = y, label = p_label),
        inherit.aes = FALSE,
        size = 5
      ) +
      facet_wrap(~ Setting, nrow = 2, ncol = 3) +
      scale_y_continuous(
        limits = c(0, 10800),
        breaks = seq(0, 10000, by = 2500),
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        x = NULL,
        y = "Split-selection frequency",
        fill = "Variable"
      ) +
      theme_bw() +
      theme(
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(vjust = -0.5, size = 17),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.width = unit(25, "pt"),
        legend.key.height = unit(25, "pt")
      )

    ggsave(
      filename = file.path(results_dir, "figure_S1_test_unbiasedness_LBRCtrees.pdf"),
      plot = p,
      width = 10,
      height = 7
    )

    cat("figure_S1_test_unbiasedness_LBRCtrees.pdf generated\n")
  }else if(mode == "figure_S11_sensitivity_analysis_unbiasedness"){
    dat <- read.csv(file.path(results_dir, "sensitivity_unbiasedness.csv"),
                    stringsAsFactors = FALSE)

    dat_long <- dat %>%
      mutate(
        Scenario = case_when(
          Scenario == "unbias_texpt" ~ "Texp",
          Scenario == "unbias_covd"  ~ "Covd",
          TRUE ~ Scenario
        ),
        mu_label = paste0("\u03bc = ", mu),
        Method = case_when(
          Method == "LTRCctree"  ~ "LTRC",
          Method == "LBRCctreeC" ~ "C",
          Method == "LBRCctreeF" ~ "F",
          TRUE ~ Method
        )
      ) %>%
      pivot_longer(
        cols = X1:X6,
        names_to = "Variable",
        values_to = "Frequency"
      ) %>%
      group_by(Scenario, mu, mu_label, Method) %>%
      mutate(
        Total = sum(Frequency),
        Percent = Frequency / Total,
        Label = percent(Percent, accuracy = 0.1)
      ) %>%
      ungroup() %>%
      mutate(
        Scenario = factor(Scenario, levels = c("Texp", "Covd")),
        mu_label = factor(
          mu_label,
          levels = paste0("\u03bc = ", sort(unique(mu)))
        ),
        Method = factor(Method, levels = c("LTRC", "C", "F")),
        Variable = factor(Variable, levels = paste0("X", 1:6))
      )

    pdat <- dat_long %>%
      distinct(Scenario, mu_label, Method, Total, pval) %>%
      mutate(
        y = Total * 1.06,
        p_label = fmt_p(pval)
      )

    p <- ggplot(dat_long, aes(x = Method, y = Frequency, fill = Variable)) +
      geom_col(width = 0.75, color = "white", linewidth = 0.5) +
      geom_text(
        aes(label = Label),
        position = position_stack(vjust = 0.5),
        size = 5
      ) +
      geom_text(
        data = pdat,
        aes(x = Method, y = y, label = p_label),
        inherit.aes = FALSE,
        size = 5
      ) +
      facet_grid(Scenario ~ mu_label) +
      scale_y_continuous(
        limits = c(0, 10800),
        breaks = seq(0, 10000, by = 2500),
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        x = NULL,
        y = "Split-selection frequency",
        fill = "Variable"
      ) +
      theme_bw() +
      theme(
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        strip.text = element_text(size = 17),
        axis.text.x = element_text(vjust = -0.5, size = 17),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.width = unit(25, "pt"),
        legend.key.height = unit(25, "pt")
      )

    ggsave(
      filename = file.path(results_dir, "figure_S11_sensitivity_analysis_unbiasedness.pdf"),
      plot = p,
      width = 15,
      height = 7
    )

    cat("figure_S11_sensitivity_analysis_unbiasedness.pdf generated\n")
  }
}


##### figure S12: sensitivity analysis on tree recovery and prediction #####
Sensitivity_plot <- function(results_dir){
  scenarios_all <- c("tree_texpt", "tree_covd", "nlin_texpt", "nlin_covd")
  Dist <- "WI"
  N <- 200
  C <- 20
  tau <- qweibull(0.9999, 2, 10) + 1
  rho_grid <- c(0, tau/10, tau/5, tau/2)

  out_rr <- list()
  out_l2 <- list()

  # ---------- Data Extraction Loop ----------
  for (scenario in scenarios_all) {
    # Determine the underlying model string
    model <- ifelse(grepl("tree", scenario), "tree", "nonlinear")

    for (rho in rho_grid) {
      trunc_rate <- round(rho/tau, 2)
      # print(paste("Scenario:", scenario, "- N:", N, "- trunc_rate:", trunc_rate))

      # Load intermediate results
      if (trunc_rate > 0) {
        result_location <- paste0(results_dir, "/sensitivity_analysis/", model, "/", Dist, "/N", N)
        setwd(result_location)

        if (scenario %in% c("tree_texpt", "nlin_texpt")) {
          fname <- sprintf("LBRC_DIST_%s_MODEL_%s_P%1.0f_N%1.0f_C%1.0f_T%1.2f", Dist, model, 30, N, C, trunc_rate)
        } else if (scenario %in% c("tree_covd", "nlin_covd")) {
          fname <- sprintf("LBRC_DIST_%s_MODEL_%s_P%1.0f_N%1.0f_C%1.0f_S%1.2f", Dist, model, 30, N, C, trunc_rate)
        }
        load(fname)
      } else {
        result_location <- paste0(results_dir, "/", model, "/", Dist, "/N", N)
        setwd(result_location)
        fname <- sprintf("LBRC_DIST_%s_MODEL_%s_P%1.0f_N%1.0f_C%1.0f", Dist, model, 30, N, C)
        load(fname)
      }

      # ---------- RR summary (Only for Tree scenarios) ----------
      if (model == "tree") {
        rr_ltrc <- mean(RES$RR$LTRCctree, na.rm = TRUE)
        rr_c    <- mean(RES$RR$LBRCctreeC, na.rm = TRUE)
        rr_f    <- mean(RES$RR$LBRCctreeF, na.rm = TRUE)

        out_rr[[length(out_rr) + 1]] <- data.frame(
          Scenario = scenario, N = N, rho = rho, mu = trunc_rate,
          method = c("LTRC", "LBRC-C", "LBRC-F"),
          RR = c(rr_ltrc, rr_c, rr_f)
        )
      }

      # ---------- L2 summary (only up to trunc_rate = 1, for all scenarios) ----------
      if (trunc_rate <= 1) {
        l2_ltrc_tree   <- RES$L2$LTRCctree
        l2_c_tree      <- RES$L2$LBRCctreeCC
        l2_f_tree      <- RES$L2$LBRCctreeFF
        l2_ltrc_forest <- RES$L2$LTRCcforest$mtryOPT
        l2_c_forest    <- RES$L2$LBRCcforestCC$mtryOPT
        l2_f_forest    <- RES$L2$LBRCcforestFF$mtryOPT

        out_l2[[length(out_l2) + 1]] <- bind_rows(
          data.frame(Scenario = scenario, N = N, rho = rho, mu = trunc_rate, method = "LTRC",   algorithm = "Tree (CIT)",   L2 = l2_ltrc_tree),
          data.frame(Scenario = scenario, N = N, rho = rho, mu = trunc_rate, method = "LBRC-C", algorithm = "Tree (CIT)",   L2 = l2_c_tree),
          data.frame(Scenario = scenario, N = N, rho = rho, mu = trunc_rate, method = "LBRC-F", algorithm = "Tree (CIT)",   L2 = l2_f_tree),
          data.frame(Scenario = scenario, N = N, rho = rho, mu = trunc_rate, method = "LTRC",   algorithm = "Forest (CIF)", L2 = l2_ltrc_forest),
          data.frame(Scenario = scenario, N = N, rho = rho, mu = trunc_rate, method = "LBRC-C", algorithm = "Forest (CIF)", L2 = l2_c_forest),
          data.frame(Scenario = scenario, N = N, rho = rho, mu = trunc_rate, method = "LBRC-F", algorithm = "Forest (CIF)", L2 = l2_f_forest)
        )
      }
    }
  }

  # ---------- Data Formatting & Factorizing ----------
  method_levels <- c("LTRC", "LBRC-C", "LBRC-F")
  method_colors <- c("LTRC" = "olivedrab3", "LBRC-C" = "dodgerblue", "LBRC-F" = "tomato")

  # 1. RR Data Formatting
  rr_df <- bind_rows(out_rr)
  rr_df$mu_f <- factor(rr_df$mu)
  rr_df$method <- factor(rr_df$method, levels = method_levels)
  rr_df$Row <- ifelse(grepl("texpt", rr_df$Scenario), "Texp", "Covd")
  rr_df$Row <- factor(rr_df$Row, levels = c("Texp", "Covd"))

  # 2. L2 Data Formatting & Facet Labeling
  l2_df <- bind_rows(out_l2)
  l2_df$mu_f <- factor(l2_df$mu)
  l2_df$method <- factor(l2_df$method, levels = method_levels)

  l2_df <- l2_df %>%
    mutate(
      ModelLabel = ifelse(grepl("tree", Scenario), "Tree", "Nonlinear"),
      RowLabel   = ifelse(grepl("texpt", Scenario), "Texp", "Covd"),
      AlgoLabel  = ifelse(grepl("CIT", algorithm), "CIT", "CIF"),
      facet_label = paste(ModelLabel, RowLabel, AlgoLabel, sep = ", ")
    )

  # Order the factors to enforce a 2x4 grid aligning top row to Texp and bottom row to Covd
  l2_df$facet_label <- factor(l2_df$facet_label, levels = c(
    "Tree, Texp, CIT", "Tree, Texp, CIF", "Nonlinear, Texp, CIT", "Nonlinear, Texp, CIF",
    "Tree, Covd, CIT",  "Tree, Covd, CIF",  "Nonlinear, Covd, CIT",  "Nonlinear, Covd, CIF"
  ))

  # ---------- Plotting Setup ----------
  shared_theme <- theme_bw() +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
      strip.text       = element_text(size = 14),
      axis.text.x      = element_text(hjust = 0.5, size = 13),
      axis.text.y      = element_text(size = 11),
      axis.title.y     = element_text(size = 16, vjust = 2),
      axis.title.x     = element_text(size = 15, vjust = -1),
      legend.position  = "none"
    )

  # ---------- Plot 1: Recovery Rate (1 Column, 2 Rows) ----------
  p_rr <- ggplot(rr_df, aes(x = mu_f, y = RR, group = method, color = method)) +
    geom_line(data = filter(rr_df, method == "LTRC"), aes(x = as.numeric(factor(mu_f)), y = RR), linewidth = 0.9, alpha = 0.8) +
    geom_line(data = filter(rr_df, method == "LBRC-C"), aes(x = as.numeric(factor(mu_f)) - 0.05, y = RR), linewidth = 0.9, alpha = 0.8) +
    geom_line(data = filter(rr_df, method == "LBRC-F"), aes(x = as.numeric(factor(mu_f)) + 0.05, y = RR), linewidth = 0.9, alpha = 0.8) +
    geom_point(data = filter(rr_df, method == "LTRC"), aes(x = as.numeric(factor(mu_f)), y = RR), size = 2.5, alpha = 0.8) +
    geom_point(data = filter(rr_df, method == "LBRC-C"), aes(x = as.numeric(factor(mu_f)) - 0.05, y = RR), size = 2.5, alpha = 0.8) +
    geom_point(data = filter(rr_df, method == "LBRC-F"), aes(x = as.numeric(factor(mu_f)) + 0.05, y = RR), size = 2.5, alpha = 0.8) +
    scale_x_discrete(limits = levels(rr_df$mu_f), labels = levels(rr_df$mu_f)) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(), labels = scales::number_format(accuracy = 0.1)) +
    scale_color_manual(values = method_colors) +
    facet_wrap(~ Row, ncol = 1) +
    labs(x = expression(mu), y = "Recovery rate") +
    shared_theme

  # ---------- Plot 2: Integrated L2 Distance (4 Columns, 2 Rows via facet_wrap) ----------
  p_l2 <- ggplot(l2_df, aes(x = mu_f, y = log(L2), fill = method)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.9) +
    scale_fill_manual(values = method_colors) +
    facet_wrap(~ facet_label, nrow = 2, ncol = 4, scales = "free_y") +
    scale_y_continuous(breaks = scales::pretty_breaks(), labels = scales::number_format(accuracy = 0.1)) +
    labs(x = expression(mu), y = "Log of integrated L2 distance") +
    shared_theme

  # ---------- Assembly & Shared Legend ----------
  # Generate an independent legend to be applied at the bottom
  legend_plot <- ggplot(rr_df, aes(x = mu_f, y = RR, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_color_manual(values = method_colors) +
    theme_bw() +
    theme(
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 14),
      legend.key.width = unit(2, "cm")
    )
  shared_legend <- get_legend(legend_plot)

  # Arrange RR on the left (1 unit) and L2 on the right (4 units)
  grid_plots <- ggarrange(
    p_rr, p_l2,
    nrow = 1,
    ncol = 2,
    widths = c(1, 4)
  )

  # Combine main grid with bottom legend
  p <- ggarrange(
    grid_plots,
    shared_legend,
    ncol = 1,
    nrow = 2,
    heights = c(1, 0.08)
  )

  # ---------- Save to File ----------
  setwd(results_dir)
  ggsave(filename = "figure_S12_sensitivity_analysis_prediction.pdf",
         plot = p, width = 15, height = 6, units = "in")

  cat("figure_S12_sensitivity_analysis_prediction.pdf generated\n")
}



generate_figures_tables <- function(mode, results_dir){
  if(mode == "figure_1_compare_recovery_rate"){
    RR_plot(results_dir)
  }else if(mode == "figure_2_LBRCtrees_vs_LTRCtree_ANOVA"){
    ANOVA_plot(results_dir)
  }else if(mode %in% c("figure_3_WI_LBRCforests_OOB_tuning_brier",
                       "figure_S2_WI_LBRCforests_OOB_tuning_cindex",
                       "figure_S3_WD_LBRCforests_OOB_tuning_brier",
                       "figure_S4_WD_LBRCforests_OOB_tuning_cindex",
                       "figure_S5_LgnBat_LBRCforests_OOB_tuning_brier",
                       "figure_S6_LgnBat_LBRCforests_OOB_tuning_cindex")){
    if(mode %in% c("figure_S2_WI_LBRCforests_OOB_tuning_cindex",
                   "figure_S4_WD_LBRCforests_OOB_tuning_cindex",
                   "figure_S6_LgnBat_LBRCforests_OOB_tuning_cindex")){
      OOBtune_plot(mode, tune.metric = "cindex", results_dir)
    }else{
      OOBtune_plot(mode, tune.metric = "brier", results_dir)
    }
  }else if(mode %in% c("figure_4_WI_20_compare_prediction_accuracy",
                       "figure_S7_WI_50_compare_prediction_accuracy",
                       "figure_S8_WD_20_compare_prediction_accuracy",
                       "figure_S9_WD_50_compare_prediction_accuracy",
                       "figure_S10_LgnBat_2050_compare_prediction_accuracy")){
    L2_plot(mode, results_dir)
  }else if(mode %in% c("figure_S1_test_unbiasedness_LBRCtrees",
                       "figure_S11_sensitivity_analysis_unbiasedness")){
    Unbiasedness_table(mode, results_dir)
    Unbiasedness_figure(mode, results_dir)
  }else if(mode == "figure_S12_sensitivity_analysis_prediction"){
    Sensitivity_plot(results_dir)
  }else{
    cat("select the mode")
  }
}









