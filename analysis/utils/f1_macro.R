library(caret)

extract_split_vars <- function(tree_struct, obj, env) {
  if (is.null(tree_struct)) {
    return(NULL)
  }

  var_id <- tree_struct$split$varid

  # If there is a split at this node
  if (!is.null(var_id)) {
    var_name <- names(obj$data)[var_id]
  } else {
    var_name <- NULL
  }

  # Recursively check the left and right child nodes for more variables
  left_var <- extract_split_vars(tree_struct$kids[[1]], obj, env)
  if (!is.null(left_var)) {
    env$left_vars <- c(env$left_vars, left_var)
  }

  right_var <- extract_split_vars(tree_struct$kids[[2]], obj, env)
  if (!is.null(right_var)) {
    env$right_vars <- c(env$right_vars, right_var)
  }

  return(var_name)
}

tree_depth <- function(tree_struct) {
  if (is.null(tree_struct$kids)) {
    return(1)  # Leaf node has a depth of 1
  } else {
    return(1 + max(sapply(tree_struct$kids, tree_depth)))
  }
}

f1_macro_tree <- function(obj, true_vars, tr_trmnd, newdata) {
  tree_struct <- node_party(obj)

  # Create a new environment to hold left_vars and right_vars
  env <- new.env()
  env$left_vars <- c()
  env$right_vars <- c()

  root_var <- extract_split_vars(tree_struct, obj, env)
  depth <- tree_depth(tree_struct)

  recovery <- 0
  f1_macro <- 0
  root_bool <- left_bool <- right_bool <- FALSE
  if (depth == 3) {
    if (root_var == "X1") root_bool <- TRUE
    if (length(env$left_vars) == 1) {
      if (env$left_vars == "X2") left_bool <- TRUE
    }
    if (length(env$right_vars) == 1) {
      if (env$right_vars == "X3") right_bool <- TRUE
    }

    if (root_bool == TRUE && left_bool == TRUE && right_bool == TRUE) {
      recovery <- 1

      tr_trmnd <- factor(tr_trmnd)
      predicted_labels <- factor(predict(obj, newdata = newdata, type = 'node'))
      conf_matrix <- confusionMatrix(predicted_labels, tr_trmnd)
      macro_prec <- mean(conf_matrix$byClass[,'Precision'])
      macro_rec <- mean(conf_matrix$byClass[,'Recall'])
      f1_macro <- 2 * (macro_prec * macro_rec) / (macro_prec + macro_rec)
    }
  }

  list(recovery = recovery, f1_macro = f1_macro)
}
