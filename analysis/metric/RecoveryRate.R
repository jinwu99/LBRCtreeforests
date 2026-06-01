extract_split_vars <- function(tree_struct, obj, env) {
  if (is.null(tree_struct)) {
    return(NULL)
  }

  var_id <- tree_struct$split$varid
  if (!is.null(var_id)) {
    var_name <- names(obj$data)[var_id]
  } else {
    var_name <- NULL
  }

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
    return(1)
  } else {
    return(1 + max(sapply(tree_struct$kids, tree_depth)))
  }
}

root_node_return <- function(obj) {
  tree_struct <- node_party(obj)
  env <- new.env()
  root_var <- extract_split_vars(tree_struct, obj, env)
  return(root_var)
}

recoverTree <- function(obj) {
  tree_struct <- node_party(obj)

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
    }
  }

  return(recovery)
}
