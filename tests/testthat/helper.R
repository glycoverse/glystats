exp_2groups <- function() {
  glyexp::filter_obs(test_gp_exp, group %in% c("C", "H"))
}

# Create a test dataset that satisfies Topliss ratio (n/p >= 5) with strong discrimination
exp_topliss_valid <- function() {
  # Create a dataset with 50 samples and 8 variables (ratio = 6.25)
  # Use deterministic patterns for maximum discrimination
  set.seed(42)  # For reproducibility

  # Generate sample names
  sample_names <- paste0("Sample_", 1:50)

  # Generate sample information (glyexp will add the "sample" column automatically)
  sample_info <- data.frame(
    group = rep(c("Control", "Treatment"), each = 25),
    stringsAsFactors = FALSE
  )
  rownames(sample_info) <- sample_names

  # Generate expression matrix (8 variables, 50 samples)
  expr_mat <- matrix(nrow = 8, ncol = 50)
  rownames(expr_mat) <- paste0("Feature_", 1:8)
  colnames(expr_mat) <- sample_names

  # Create highly discriminative patterns
  # Control group (columns 1-25): specific pattern
  for (i in 1:25) {
    # Features 1-4: high expression in Control
    expr_mat[1:4, i] <- c(8, 9, 7.5, 8.5) + rnorm(4, 0, 0.2)
    # Features 5-8: low expression in Control
    expr_mat[5:8, i] <- c(2, 1.5, 2.5, 1.8) + rnorm(4, 0, 0.2)
  }

  # Treatment group (columns 26-50): opposite pattern
  for (i in 26:50) {
    # Features 1-4: low expression in Treatment
    expr_mat[1:4, i] <- c(2, 1.5, 2.5, 1.8) + rnorm(4, 0, 0.2)
    # Features 5-8: high expression in Treatment
    expr_mat[5:8, i] <- c(8, 9, 7.5, 8.5) + rnorm(4, 0, 0.2)
  }

  # Ensure all values are positive (for log transformation)
  expr_mat <- pmax(expr_mat, 0.1)

  # Create variable information
  var_info <- data.frame(
    variable_type = rep("biomarker", 8),
    stringsAsFactors = FALSE
  )
  rownames(var_info) <- rownames(expr_mat)

  # Create glyexp object
  glyexp::experiment(expr_mat = expr_mat, sample_info = sample_info, var_info = var_info,
                     exp_type = "glycomics", glycan_type = "N")
}

# Create a multi-group dataset that satisfies Topliss ratio
exp_multigroup_valid <- function() {
  # Create a dataset with 40 samples and 6 variables (ratio = 6.67)
  set.seed(456)  # For reproducibility

  # Generate sample information with 4 groups (glyexp will add the "sample" column automatically)
  sample_info <- data.frame(
    group = rep(c("A", "B", "C", "D"), each = 10),
    stringsAsFactors = FALSE
  )
  rownames(sample_info) <- paste0("S", 1:40)

  # Generate expression matrix (6 variables, 40 samples)
  expr_mat <- matrix(nrow = 6, ncol = 40)
  rownames(expr_mat) <- paste0("Gene", 1:6)
  colnames(expr_mat) <- rownames(sample_info)

  # Create group-specific patterns
  for (i in 1:4) {
    start_col <- (i-1)*10 + 1
    end_col <- i*10

    # Each group has different expression patterns
    expr_mat[1:2, start_col:end_col] <- matrix(rnorm(2*10, mean = i+2, sd = 0.8), nrow = 2)
    expr_mat[3:4, start_col:end_col] <- matrix(rnorm(2*10, mean = 6-i, sd = 0.8), nrow = 2)
    expr_mat[5:6, start_col:end_col] <- matrix(rnorm(2*10, mean = 3, sd = 0.5), nrow = 2)
  }

  # Ensure all values are positive
  expr_mat <- pmax(expr_mat, 0.1)

  # Create variable information
  var_info <- data.frame(
    variable_type = rep("gene", 6),
    stringsAsFactors = FALSE
  )
  rownames(var_info) <- rownames(expr_mat)

  # Create glyexp object
  glyexp::experiment(expr_mat = expr_mat, sample_info = sample_info, var_info = var_info,
                     exp_type = "glycomics", glycan_type = "N")
}
