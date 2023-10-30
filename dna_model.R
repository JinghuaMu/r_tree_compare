# +
dna_model <- function(model_string) {
  # Store two pairs of vectors
  model_to_code <- list(
  "JC" = "000000",
  "F81" = "000000",
  "K80" = "010010",
  "K2P" = "010010",
  "HKY" = "010010",
  "HKY85" = "010010",
  "TN" = "010020",
  "TN93" = "010020",
  "TNe" = "010020",
  "K81" = "012210",
  "K3P" = "012210",
  "K3Pu" = "012210",
  "K81u" = "012210",
  "TPM2" = "010212",
  "TPM2u" = "010212",
  "TPM3" = "012012",
  "TPM3u" = "012012",
  "TIM" = "012230",
  "TIMe" = "012230",
  "TIM2" = "010232",
  "TIM2e" = "010232",
  "TIM3" = "012032",
  "TIM3e" = "012032",
  "TVM" = "012314",
  "TVMe" = "012314",
  'SYM' = '012345', 
 'GTR'= '012345'
 )

model_to_if_equal_base_freq <- list(
   'JC' = TRUE,
   'F81' = FALSE,
   'K80' = TRUE,
   'K2P' = TRUE,
   'HKY' = FALSE,
   'HKY85' = FALSE,
   'TN' = FALSE,
   'TN93' = FALSE,
   'TNe' = TRUE,
   'K81' = TRUE,
   'K3P' = TRUE,
   'K3Pu' = FALSE,
   'K81u' = FALSE,
   'TPM2' = TRUE,
   'TPM2u' = FALSE,
   'TPM3' = TRUE,
   'TPM3u' = FALSE,
   'TIM' = FALSE,
   'TIMe' = TRUE,
   'TIM2' = FALSE,
   'TIM2e' = TRUE,
   'TIM3' = FALSE,
   'TIM3e'=TRUE,
   'TVM'=FALSE, 
    'TVMe'=TRUE, 
    'SYM'=TRUE, 
    'GTR'=FALSE
 )

 # Split the input string by ";"
 models <- strsplit(model_string, ";")[[1]]
 results <- tibble()

 for (model_string in models) {
   
    # Split the input string by "+"
    split_string <- strsplit(model_string, "\\+")
    
    # Get the model and parameters
    model_params <- strsplit(split_string[[1]][1], "\\{")
    
    model <- gsub("\\{.*", "", model_params[[1]][1])
    
    params <- if (length(model_params[[1]]) > 1) {
      as.numeric(strsplit(gsub("\\}", "", model_params[[1]][2]), ",")[[1]])
    } else {
      NA
    }
    
    # Get the frequence type and vaues
    F_type <- gsub("(.*)\\{.*", "\\1", split_string[[1]][2])
    if (is.na(F_type)) {
        F_type <- "FQ"
    }
    base_freqs <- if (length(split_string[[1]]) > 1) {
      if (grepl("\\}", split_string[[1]][2])) {
          as.numeric(strsplit(gsub(".*\\{([^}]*)\\}.*", "\\1", split_string[[1]][2]), ",")[[1]])
      } else {
        rep(0.25, 4)  
      }
    } else {
      rep(0.25, 4)
    }

    # Check the code and if_equal_base_freq
    code <- model_to_code[[model]]
    if_equal_base_freq <- model_to_if_equal_base_freq[[model]]
    
     # Calculate R
     R <- if (!any(is.na(params))) {
       if(length(params) == 6){
         params
       } else {
         params <- c(1, params)
         sapply(strsplit(code, "")[[1]], function(x) params[as.integer(x) + 1])
       }
     } else {
       params <- c(1)
       sapply(strsplit(code, "")[[1]], function(x) params[as.integer(x) + 1])
     }
     
     # Calculate F
     F <- if (if_equal_base_freq) {
       rep(1/4, 4)
     } else if (!any(is.na(base_freqs))) {
       base_freqs
     } else {
       rep(NA, 4)
     }
    
    names(R) <- c("A-C", "A-G", "A-T", "C-G", "C-T", "G-T")
    names(F) <- c("A", "C", "G", "T")
    
    # Function to get rate
    get_rate <- function(i, j, R) {
      if (paste(names[i], names[j], sep = "-") %in% names(R)) {
        return(R[paste(names[i], names[j], sep = "-")])
      } else {
        return(R[paste(names[j], names[i], sep = "-")])
      }
    }

    # Calculate Q
    Q <- matrix(0, nrow = 4, ncol = 4)
    names <- c("A", "C", "G", "T")
    rownames(Q) <- names
    colnames(Q) <- names

    for (i in 1:4) {
      for (j in 1:4) {
        if (i != j) {
          Q[i, j] <- F[names[j]] * get_rate(i, j, R)
          }
      }
      Q[i, i] <- -sum(Q[i, ])
    }

      # Rescale Q so that the mean rate of substitution is one
    mu <- -1 / sum(F * diag(Q))
    Q <- mu * Q
    
    # Add to results
    results <- bind_rows(results, tibble(
      class_id = nrow(results) + 1,
      class = model,
      F_type = F_type,
      df_R = length(unique(strsplit(code, "")[[1]])) - 1,
      R = list(R),
      F = list(F),
      Q = list(Q),
      base_freq = list(base_freqs)
    ))
  }
  return(results)
}

