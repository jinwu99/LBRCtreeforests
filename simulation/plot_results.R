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
          panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
          strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
          strip.text       = element_text(size = 14),
          axis.text.x      = element_text(hjust = 0.5, size=11),
          axis.text.y      = element_text(size=12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          legend.text     = element_text(size = 12),

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

  all_dists_with_legend <- ggarrange(
    all_dists,
    get_legend(legend_plot),
    ncol = 1,
    heights = c(1, 0.05)
  )

  p <- annotate_figure(
    all_dists_with_legend,
    left = text_grob("Recovery Rate", rot = 90, size = 15)
  )

  setwd(results_dir)
  ggsave(filename = "figure_1_compare_recovery_rate.pdf",
         plot = p, width = 12, height = 6, units = "in")
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
            result_location <- paste0(results_dir,"/",model,"/",dist,"/N",N)
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
      mutate(diff = (baseline_L2-L2) )

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
        strip.text       = element_text(size = 14),
        axis.text.x      = element_text(hjust = 1, size=13),
        axis.text.y      = element_text(size = 12),
        panel.spacing.x  = unit(0.3, "lines"),
        axis.title.y = element_text(size = 13)
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
        size      = 3.5,
        inherit.aes = FALSE
      )

    return(gg)
  }

  ggtree <- anova_model_plot("tree")
  gglinear <- anova_model_plot("linear")
  ggnonlinear <- anova_model_plot("nonlinear")

  p <- ggtree / gglinear / ggnonlinear +
    plot_layout(heights = c(1, 1, 1))

  setwd(results_dir)
  ggsave(filename = paste0(mode,".pdf"),
         plot = p, width = 12, height = 8, units = "in")
}


##### figure 3: validation of OOB tuning for LBRCforests #####
OOBtune_plot <- function(mode, results_dir){
  if(mode == "figure_S1_LBRCforests_OOB_tuning"){
    dists <- c("Lgn","Bat")
    models <- c("tree")
    width <- 10; height <- 6
    nrow <- 2; ncol <- 3
    desired_order <- c(
      "Lgn, tree, n = 100",
      "Lgn, tree, n = 200",
      "Lgn, tree, n = 400",
      "Bat, tree, n = 100",
      "Bat, tree, n = 200",
      "Bat, tree, n = 400"
    )
  }else{ # "figure_3_LBRCforests_OOB_tuning"
    dists <- c("WI", "WD")
    models <- c("tree","linear","nonlinear")
    width <- 17; height <- 8
    nrow <- 3; ncol <- 6
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

  th <- 500

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
        tmp$model <- model
        tmp$dist <- dist
        tmp$N <- N

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

  mtrys <- c(paste(mtry_candidates),"Tuned")

  plot_df <- storage %>%
    left_join(method_lookup, by = "method") %>%
    mutate(
      mtry = factor(mtry, levels = mtrys),
      facet_label = paste0(dist, ", ", model, ", n = ", N)
    )

  plot_df <- plot_df %>%
    mutate(facet_label = factor(facet_label, levels = desired_order))

  p <- ggplot(plot_df, aes(x = mtry, y = mean, fill = Model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    facet_wrap(~ facet_label,
               nrow = nrow,
               ncol = ncol,
               scales = "free_y") +
    scale_fill_manual(
      values = c(
        "LTRC-CIF"   = "olivedrab3",
        "LBRC-CIF-C"  = "dodgerblue",
        "LBRC-CIF-F"   = "tomato"
      )
    ) +
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
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    )

  setwd(results_dir)
  ggsave(filename = paste0(mode,".pdf"),
         plot = p, width = width, height = height, units = "in")
}


##### figure 4: prediction comparison across models #####
L2_plot <- function(mode, results_dir){
  if(mode == "figure_S2_compare_prediction_accuracy"){
    dists <- "WD"
    models <- c("tree","linear","nonlinear")
    width <- 17; height <- 8
  }else if(mode == "figure_S3_compare_prediction_accuracy"){
    dists <- c("Lgn", "Bat")
    models <- "tree"
    width <- 17; height <- 6
  }else{ # "figure_4_compare_prediction_accuracy"
    dists <- "WI"
    models <- c("tree","linear","nonlinear")
    width <- 17; height <- 8
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
        for(C in c(20,50)) {
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
    facet_wrap(~ facet_label, nrow = 3, ncol=6, scales = "free_y") +
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
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    )

  setwd(results_dir)
  ggsave(filename = paste0(mode,".pdf"),
         plot = p, width = width, height = height, units = "in")
}


plot_figures <- function(mode, results_dir){
  if(mode == "figure_1_compare_recovery_rate"){
    RR_plot(results_dir)
  }else if(mode == "figure_2_LBRCtrees_vs_LTRCtree_ANOVA"){
    ANOVA_plot(results_dir)
  }else if(mode %in% c("figure_3_LBRCforests_OOB_tuning",
                       "figure_S1_LBRCforests_OOB_tuning")){
    OOBtune_plot(mode, results_dir)

  }else if(mode %in% c("figure_4_compare_prediction_accuracy",
                       "figure_S2_compare_prediction_accuracy",
                       "figure_S3_compare_prediction_accuracy")){
    L2_plot(mode, results_dir)
  }
}









