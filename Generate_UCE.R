# function for generating a sequence
generate_sequence <- function(length) {
  return(paste(sample(c("A", "T", "G", "C"), length, replace = TRUE), collapse = ""))
}

# function for mutating a sequence
mutate_sequence <- function(sequence, mutation_rate) {
  sequence <- strsplit(sequence, "")[[1]]
  mutate_positions <- sample(1:length(sequence), mutation_rate * length(sequence))
  for (position in mutate_positions) {
    sequence[position] <- sample(c("A", "T", "G", "C"), 1)
  }
  return(paste(sequence, collapse = ""))
}

# function for padding a sequence
pad_sequence <- function(sequence, length, dir = "L") {
    if (dir == "L") {
       return(paste(paste(rep("-", length), collapse = ""),sequence, sep = "")) 
    }else {
       return(paste(sequence, paste(rep("-", length), collapse = ""), sep = ""))
    }
}

# function for introducing directed mutations at a specific position
directed_mutate_sequence <- function(sequences, position, target_base, mutation_base, mutation_rate) {
  # iterate over each sequence
  for (i in 1:length(sequences)) {
    sequence <- strsplit(sequences[i], "")[[1]]
    # check if the base at the position is the target base
    if (sequence[position] == target_base) {
      # introduce mutation with the specified rate
      if (runif(1) < mutation_rate) {
        sequence[position] <- mutation_base
      }
    }
    # join the sequence back together
    sequences[i] <- paste(sequence, collapse = "")
  }
  return(sequences)
}


# generate core and flanks
core <- "ATGAAAAGGCTTGAGTGAAG"
left_flank <- generate_sequence(60)
right_flank <- generate_sequence(60)

# duplicate the regoin of core
sequences <- replicate(15, core)

# +
strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
# Add flanks for each sequence
for (i in 1:length(sequences)) {
  left_length <- pmax(1, pmin(round(rnorm(1, mean = 10, sd = 0)), 60))
  right_length <- pmax(1, pmin(round(rnorm(2, mean = 10, sd = 0)), 60))

  left_sequence <- strReverse(substring(left_flank, 1, left_length))
  right_sequence <- substring(right_flank, 1, right_length)

  left_sequence <- mutate_sequence(left_sequence, mutation_rate = 0.4)
  right_sequence <- mutate_sequence(right_sequence, mutation_rate = 0.4)

  sequences[i] <- paste(left_sequence, 
                        mutate_sequence(sequences[i], mutation_rate = 0.05), 
                        right_sequence, sep = "")
}

sequences <- directed_mutate_sequence(sequences, 20, "C", "A", 0.5)
sequences <- directed_mutate_sequence(sequences, 20, "G", "A", 0.5)
sequences <- directed_mutate_sequence(sequences, 20, "T", "A", 0.5)
# -

# write the sequences to fasta
fasta_file <- file("/home/tim/mixmodel/R_tree_comparison/sequences_region.fasta", "w")
for (i in 1:length(sequences)) {
  cat(sprintf(">seq%d\n", i), file = fasta_file)
  cat(sequences[i], "\n", file = fasta_file)
}
close(fasta_file)




