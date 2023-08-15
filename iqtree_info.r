summarise_iqtree <- function(file_path) {
  # open iqtree file
  iqtree_text <- readLines(file_path)
  # abstract iqtree info
  iqtree_summarise <- list()

  # the best model and its criteria
  line_best_model <- grep("Best-fit model according to", iqtree_text)
  iqtree_summarise$model_crireia <- sub(".*to\\s+(.*):.*", "\\1", iqtree_text[line_best_model])
  iqtree_summarise$best_model <- sub(".*:\\s*(.*)", "\\1", iqtree_text[line_best_model])
  line_outgroup <- grep("Tree is UNROOTED", iqtree_text)
  iqtree_summarise$outgroup <- sub(".*'\\s*(.*)'.*", "\\1", iqtree_text[line_outgroup])
  rm(line_best_model, line_outgroup)

  # Find the line numbers for the substitution model and likelihood tree
  line_sub <- grep("SUBSTITUTION PROCESS", iqtree_text)
  line_likelihood <- grep("MAXIMUM LIKELIHOOD TREE", iqtree_text)

  # Extract the substitution model information
  model_info <- iqtree_text[(line_sub + 3):(line_likelihood - 1)]
  model_substitution <- sub(".*: (.*)", "\\1", model_info[1])
  model_components <- unlist(strsplit(model_substitution, "\\+"))

  # Check if the model includes a rate heterogeneity mixture
  if (grepl("(m|M)ixture", model_info[1])) {
    # Find the line numbers for the rate heterogeneity and site proportion and rates sections
    line_RHAS <- grep("Model of rate heterogeneity:", model_info)

    # Extract the rate heterogeneity and site proportion and rates information
    iqtree_summarise$model_info <- read.table(text = model_info[3:(line_RHAS - 1)], header = TRUE)
    if (grepl("R", model_components[length(model_components)])) {
      block_SPR <- model_info[(line_RHAS + 2):(length(model_info) - 1)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      rm(block_SPR)
    } else if (grepl("G", model_components[length(model_components)])) {
      block_SPR <- model_info[(line_RHAS + 2):(length(model_info) - 2)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      rm(block_SPR)
    } else {
      iqtree_summarise$RHAS <- NULL
    }
  } else {
    # Find the line numbers for the rate parameter R, state frequencies, rate matrix Q, rate heterogeneity, and site proportion and rates sections
    line_R <- grep("Rate parameter R", model_info)
    line_F <- grep("State frequencies", model_info)
    line_Q <- grep("Rate matrix Q", model_info)
    line_RHAS <- grep("Model of rate heterogeneity:", model_info)

    # Extract the rate parameter R, state frequencies, and rate matrix Q information
    block_R <- model_info[(line_R + 2):(line_F - 2)]
    block_Q <- model_info[(line_Q + 2):(line_RHAS - 1)]

    iqtree_summarise$R <- as.numeric(sub(".*: ([0-9.]+)", "\\1", block_R))
    names(iqtree_summarise$R) <- sub("(.*): [0-9.]+", "\\1", block_R)

    if (grepl("equal frequencies", model_info[line_F])) {
      iqtree_summarise$F <- rep(1 / 4, 4)
    } else {
      block_F <- model_info[(line_F + 2):(line_Q - 2)]
      iqtree_summarise$F <- as.numeric(sub(".*= ([0-9.]+)", "\\1", block_F))
      rm(block_F)
    }
    names(iqtree_summarise$F) <- c("A", "C", "G", "T")

    iqtree_summarise$Q <- read.table(text = block_Q, header = FALSE, row.names = 1)
    colnames(iqtree_summarise$Q) <- rownames(iqtree_summarise$Q)

    # Create a parameter string for the model
    parameter <- paste0(
      model_components[1], "{", paste(as.character(iqtree_summarise$R), collapse = ","), "}+",
      model_components[2], "{", paste(as.character(iqtree_summarise$F), collapse = ","), "}"
    )

    # Store the model information in a data frame
    iqtree_summarise$model_info <- data.frame(
      "No" = 1,
      "Component" = model_components[1],
      "Rate" = 1.0000,
      "Weight" = 1.0000,
      "Parameters" = parameter
    )


    # Extract the rate heterogeneity information if present
    if (grepl("R", model_components[length(model_components)])) {
      block_SPR <- model_info[(line_RHAS + 2):(length(model_info) - 1)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      rm(block_SPR)
    } else if (grepl("G", model_components[length(model_components)])) {
      block_SPR <- model_info[(line_RHAS + 2):(length(model_info) - 2)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      rm(block_SPR)
    } else {
      iqtree_summarise$RHAS <- NULL
    }

    # Remove unnecessary variables
    rm(line_sub, line_R, line_F, line_Q, line_RHAS)
    rm(block_R, block_Q)
  }

  # tree length
  line_length <- grep("Total tree length", iqtree_text)
  iqtree_summarise$total_tree_length <- sub(".*: ([0-9.]+)", "\\1", iqtree_text[line_length])
  iqtree_summarise$internal_tree_length <- sub(".*: ([0-9.]+).*", "\\1", iqtree_text[line_length + 1])
  rm(line_length)

  # statistics

  iqtree_summarise$log_likelihood <- sub(".*: ([0-9.-]+) .*", "\\1", iqtree_text[line_likelihood + 3])
  iqtree_summarise$log_likelihood_sd <- sub(".*: ([0-9.-]+) \\(s.e. ([0-9.]+)\\)", "\\2", iqtree_text[line_likelihood + 3])
  iqtree_summarise$log_likelihood_Unconstrained <- sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 4])
  iqtree_summarise$num_free_params <- sub(".*: ([0-9]+)", "\\1", iqtree_text[line_likelihood + 5])
  iqtree_summarise$AIC <- sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 6])
  iqtree_summarise$AICc <- sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 7])
  iqtree_summarise$BIC <- sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 8])
  # iqtree_summarise$num_params <- sub(".*: ([0-9]+)", "\\1", iqtree_text[line_likelihood + 11])
  # iqtree_summarise$samp_size <- sub(".*: ([0-9]+)", "\\1", iqtree_text[line_likelihood + 12])
  rm(line_likelihood)

  return(iqtree_summarise)
}

short_description <- function(summarise_list) {
  cat("Best model: ", summarise_list$best_model, "\n")
  cat("Tree length: ", summarise_list$total_tree_length, "(inner: ", summarise_list$internal_tree_length, ")\n", sep = "")
  cat("Log-likelihood: ", summarise_list$log_likelihood, "\n")
  cat("Number of free parameters: ", summarise_list$num_free_params, "\n")
  cat("AIC, AICc, BIC: ", summarise_list$AIC, summarise_list$AICc, summarise_list$BIC, "\n")
  cat("Model parameters: \n")
  print(summarise_list$model_info$Parameters)
  cat("Rate heterogeneity: \n")
  if (nrow(summarise_list$RHAS) > 0) {
    cat("    Rate:", summarise_list$RHAS$Relative_rate, "\n")
    cat("    Prop:", summarise_list$RHAS$Proportion, "\n")
  } else {
    cat("    None\n")
  }
}
