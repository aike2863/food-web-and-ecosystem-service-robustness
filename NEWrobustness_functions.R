## functions to run robustness analysis

# functions rely on tidyverse, "purrr" in particular
library(tidyverse)

# the recursive species extinction function --------
# depends on itself and two additional functions: whos_next and decendr
# it requires three inputs: the named square interaction matrix N, the degree order
# of species removal ("h" = highest first, "l" = lowest first, "r" = randomly,
# or "g" = given order from the vector order),
# and a vector of species that will be removed (the "basal" species)
collapse <- function (N, deg, basal) {
  if (ncol(N) < 2) {
    return(list(matrix = N))
  }
  i <- whos_next(N, deg, basal)
  N2 <- descendr(N, i)
  
  recurse <- collapse(N2, deg, basal)
  
  output <- c(list(i = i, desc = colnames(N)[colnames(N) %ni% colnames(N2)], matrix = N2), recurse)
  output
}

# function to determine order of species removal ------------
# this function is written to do the classic robustness analysis where species are
# removed in order of degree, and this function determines who is removed next in
# the sequence. i added the ability to tell it the order of species removal

# whos_next <- function(N, deg, basal) {
#   degr <- rowSums(N)
#   degr_vec <- degr[names(degr) %in% basal]
#   
#   if (deg == "h") {
#     return(names(degr_vec)[order(degr_vec, 
#                                  decreasing = TRUE)][1])
#   } else if (deg == "l") {
#     return(rev(names(degr_vec)[order(degr_vec, 
#                                      decreasing = TRUE)])[1])
#   } else {
#     ifelse(length(degr_vec) < 1, return(NA),
#            return(sample(x = names(degr_vec), size = 1)))
#   }
# }

whos_next <- function(N, deg, basal) {
  degr <- rowSums(N)
  degr_vec <- degr[names(degr) %in% basal]
  
  if (deg == "h") {
    return(names(degr_vec)[order(degr_vec, 
                                 decreasing = TRUE)][1])
  } else if (deg == "l") {
    return(rev(names(degr_vec)[order(degr_vec, 
                                     decreasing = TRUE)])[1])
  } else if (deg == "g") {
    
    # this sucks up time because recalculate this each run,
    # could improve by pre-allocating
    col_N <- colnames(N)
    basal_order <- data.frame(b_order = seq(1:length(basal)),
                              basal = basal,
                              stringsAsFactors = FALSE)
    basal_filter <- filter(basal_order, basal %in% col_N)
    basal_return <- which(basal_filter$b_order == min(basal_filter$b_order))

    return(basal_filter$basal[basal_return])
    
  } else {
    ifelse(length(degr_vec) < 1, return(NA),
           return(sample(x = names(degr_vec), size = 1)))
  }
}

# function to determine who goes extinct when one species removed ----------
# pulls all the "descendent" species from the matrix - removes only those
# species who have no other resources after the focal species is removed
descendr <- function(N, i) {
  
  if (length(i) < 1) {
    return(N)
  }
  
  if (length(i) == ncol(N) |
      length(i) == ncol(N)-1) {
    return(matrix(0, nrow=1, ncol=1))
  }
  
  N2 <- N[, -c(which(colnames(N) %in% i))]
  N2 <- N2[-c(which(rownames(N2) %in% i)), ]
  i <- names(which(colSums(N2) < 1))[names(which(colSums(N2) < 1)) %ni% names(which(colSums(N) < 1))]
  
  descendr(N = N2, i = i)
  
}

# wrangle output from 'collapse' --------
# the output from the "collapse" function is long and unwieldy, it basically stores 
# everything that happens at every removal step, but generally you don't need to look
# at that. this function acts as a wrapper to the output, and works in 
# conjunction with the collapse_wrap function below
collapse_output_with_names <- function(output, basal) {
  
  # pull out all the matrices
  output_mat <- output[names(output) == "matrix"]
  # pull out all the lost taxa ("desc")
  output_desc <- output[names(output) == "desc"]
  
  # make sure all the empty rows/columns are pulled from the matrices
  output_mat <- map(output_mat, rm_empty)
  
  # make sure nothing is NULL because that causes problems
  output_i <- output[names(output) == "i"]
  output_i[sapply(output_i, is.null)] <- NA
  
  df <- data.frame(basal_removed_each = rep(NA, 1 + length(output_i)),
                   # what were the basal taxa removed at each step
                   num_lost_each = c(0, lengths(output_desc)),
                   # how many taxa were lost as the result of removing said basal taxa
                   stringsAsFactors = FALSE)
  df$basal_removed_each <- c(NA, unlist(output_i))
  
  # how many taxa remain at each step - calculate from the matrix size
  df$num_remain_tot <- c(
    # for the first, sum the remaining taxa + the removed 
    unlist(lapply(output_mat[1], ncol)) + lengths(output_desc)[1],
    # must remove the last entry from the collapse function
    # it returns 1 extra empty matrix
    unlist(lapply(output_mat[-length(output_mat)], ncol)))
  
  # how many basal taxa have been cumulatively removed
  df$num_removed_tot <- c(0:(nrow(df)-1))
  # how many taxa have been cumulatively lost
  df$num_lost_tot <- cumsum(df$num_lost_each)
  
  # calculate proportions
  df$prop_remain <- df$num_remain_tot / df$num_remain_tot[1]
  df$prop_removed <- df$num_removed_tot / length(basal)
  
  # add a column with names of species lost at each step
  output_names <- unlist(lapply(output[names(output) == "desc"], FUN = function(x) paste(x, collapse = ";")))
  df$name_lost <- c("", as.character(output_names))
  
  return(df)
}

# automate running a series of 'collapse' simulations --------
# this is a wrapping function that runs the 3 simulated degree removal types and 
# returns a dataframe that can directly be plotted
collapse_wrap <- function(N = N, basal = basal, n_rand = 99) {
  raw_list <- vector(mode = "list", n_rand)
  for (i in 1:n_rand) {
    raw_list[[i]] <- collapse_output_with_names(collapse(N = N, deg = c("r"), 
                                              basal = basal), basal = basal)
  }
  names(raw_list) <- paste("rand", 1:n_rand, sep = "_")
  raw_list$high <- collapse_output_with_names(collapse(N = N, deg = c("h"), basal = basal), basal = basal)
  raw_list$low <- collapse_output_with_names(collapse(N = N, deg = c("l"), basal = basal), basal = basal)
  
  raw_list %>%
    bind_rows(.id = "id") %>%
    separate(id, into = c("degree_type", "rep"), sep = "_", remove = FALSE) %>%
    return()
}

# a version of the above function that runs without any randomizations
collapse_wrap_nr <- function(N = N, basal = basal) {
  raw_list <- vector(mode = "list", 2)
  names(raw_list) <- c("high", "low")
  raw_list$high <- collapse_output_with_names(collapse(N = N, deg = c("h"), basal = basal), basal = basal)
  raw_list$low <- collapse_output_with_names(collapse(N = N, deg = c("l"), basal = basal), basal = basal)
  
  raw_list %>%
    bind_rows(.id = "degree_type") %>%
    return()
}

# run with only the "high" degree species removed first
collapse_wrap_high <- function(N = N, basal = basal) {
  tmp <- collapse_output_with_names(collapse(N = N, deg = c("h"), basal = basal), basal = basal)
  
  tmp %>%
    add_column(degree_type = "high") %>%
    return()
}

# run with species removed in a specific order
collapse_wrap_given <- function(N = N, basal = basal) {
  tmp <- collapse_output_with_names(collapse(N = N, deg = c("g"), basal = basal), basal = basal)
  
  tmp %>%
    add_column(degree_type = "given") %>%
    return()
}

# calculate area under the curve --------
robust_auc <- function(x, y) sum(diff(x) * (head(y,-1)+tail(y,-1)))/2

# calculate auc on the output from "collapse_wrap", when more than one
# simulation run
auc_wrapper <- function(output) {
  output %>%
    unite(id, degree_type, rep, col = "unique_id", sep = "-") %>%
    group_by(unique_id) %>%
    summarise(auc = robust_auc(prop_removed, prop_remain)) %>%
    separate(col = unique_id, into = c("id", "degree_type", "rep"), sep = "-") %>%
    mutate(rep = as.numeric(rep)) %>% 
    return()
}

# write new auc wrapper function for given sequence
auc_wrapper_given <- function(output_given) {
  output_given %>%
    summarise(auc = robust_auc(prop_removed, prop_remain)) %>%
    return()
}

# extra functions internal to the robustness functions --------

'%ni%' <- Negate('%in%')

# remove rows and columns that are both 0 (= 'degenerate' species)
rm_empty <- function (x) {
  tmp_x <- abs(x)
  tmp_x[is.na(tmp_x)] <- 0
  nums <- which(colSums(tmp_x) == 0 & rowSums(tmp_x) == 0)
  if(length(nums) < 1) {
    return (x)} else {return(x[-c(nums), -c(nums)])}
}
