# produce the simulation data
simulation <- function() {
  PATHWAY <- 5; GENE_NUMBER <- 15;
  
  HIGH_RELATION <- c(0.5, 0.8); LOW_RELATION <- c(0.1, 0.3);
  
  SUB_GENES <- c(5:15); SUB_GOURPS_PATHWAY <- 3;
  SUB_GOURPS_MAP <- list(c(1, 1), c(3, 2), c(4, 3));
  EXTRA_HIGH_RELATION <- c(1, 2, 3, 4);
  
  
  RANK <- 500; MEAN <- rep(0, RANK);
  DATASET_SIZE <- 1000;
  
  CONSTRUCT_NOISE_INTERVAL <- c(0, 0.1);
  
  NOISE_STD <- 0.1;
  
  
  is_out_bound <- function(vector, index) {
    return(index > length(vector));
  };
  
  divide_groups <- function(genes, subgroup_number) {
    result <- list();
    genes_length <- length(genes);
    
    for (gene_index in 1: genes_length) {
      group_index <- (gene_index - 1) %% subgroup_number + 1;
      if (is_out_bound(result, group_index)) {
        result[[group_index]] <- genes[gene_index];
      } else {
        result[[group_index]] <- c(
          result[[group_index]], genes[gene_index]);
      };
    };
    return(result);
  };
  
  get_relation_groups <- function(
    sub_genes=SUB_GENES,
    sub_groups_number=SUB_GOURPS_PATHWAY,
    mapping=SUB_GOURPS_MAP,
    extra=EXTRA_HIGH_RELATION) {
    
    max_index <- 0;
    for (map in mapping) {
      max_index <- max(c(max_index), map[[1]]);
    };
    strong_relation_groups <- list();
    for (index in 1: max_index) {
      strong_relation_groups[[index]] <- c(index);
    };
    rest_groups <- divide_groups(sub_genes, sub_groups_number);
    for (mapping_item in mapping) {
      strong_relation_groups[[mapping_item[1]]] <-
        c(
          strong_relation_groups[[mapping_item[1]]],
          rest_groups[[mapping_item[2]]]
        );
    };
    strong_relation_groups[[max_index + 1]] <- extra;
    
    return(strong_relation_groups);
  };
  
  is_strong_relation_group <- function(gene1, gene2) {
    strong_relation_groups <- get_relation_groups();
    for (strong_relation_group in strong_relation_groups) {
      if ((gene1 %in% strong_relation_group)
          &&
          (gene2 %in% strong_relation_group)) {
        return(TRUE);
      };
    };
    return(FALSE);
  };
  
  construct_a_relation_matrix <- function(
    nrow=GENE_NUMBER, ncol=GENE_NUMBER,
    strong_relation_cof=HIGH_RELATION,
    low_relation_cof=LOW_RELATION) {
    
    result <- matrix(data = 0, nrow = nrow, ncol = ncol);
    for (row_index in 1: nrow) {
      for (col_index in row_index: ncol) {
        if (row_index == col_index) {
          result[row_index, col_index] <- 1;
          next;
        };
        if (is_strong_relation_group(row_index, col_index)) {
          result[row_index, col_index] <-
            runif(1, strong_relation_cof[1], strong_relation_cof[2]);
          result[col_index, row_index] <-
            result[row_index, col_index];
        } else {
          result[row_index, col_index] <-
            runif(1, low_relation_cof[1], low_relation_cof[2]);
          result[col_index, row_index] <-
            result[row_index, col_index];
        };
      };
    };
    return(result);
  };
  
  random_matrix <- function(row, col, interval) {
    matrix_data <- matrix(data = 0, nrow = row, ncol = col);
    for (row_index in 1: row) {
      for (col_index in row_index: col) {
        col_vector <- runif(1, interval[1], interval[2]);
        matrix_data[row_index, col_index] <- col_vector;
        matrix_data[col_index, row_index] <- col_vector;
      };
    };
    return(matrix_data);
  };
  
  construct_relation_matrix <- function(
    row_number=RANK, col_number=RANK,
    child_rank=GENE_NUMBER,
    construct_noise_interval=CONSTRUCT_NOISE_INTERVAL,
    number=PATHWAY) {
    
    matrix_data <- random_matrix(
      row_number, col_number, construct_noise_interval);
    diag(matrix_data) <-
      diag(diag(1, nrow = row_number, ncol = col_number));
    
    for (row_index in seq(1, child_rank * number, child_rank)) {
      matrix_tmp <- construct_a_relation_matrix();
      end_index <- (row_index + child_rank - 1);
      matrix_data[row_index: end_index, row_index: end_index] <-
        matrix_tmp;
    };
    
    # The product of a symmetric matrix must be a positive definite matrix
    return(matrix_data %*% t(matrix_data));
  };
  
  simulation_data <- function(
    sigma=construct_relation_matrix(),
    mean=MEAN, number=DATASET_SIZE,
    noise_std=NOISE_STD) {
    data <- mvrnorm(number, mean, sigma);
    
    data_cols_number <- length(data) / nrow(data);
    for (col_index in 1: data_cols_number) {
      noise <- rnorm(number, 0, noise_std);
      data[, col_index] <- data[, col_index] + noise;
    };
    return(t(data));
  };
  
  data <- simulation_data();
  return(data);
}

# cut seed genes from genes
cut_seed_genes <- function(genes, seed_gene_indexs) {
  seed_genes <- matrix(
    data = 0,
    nrow = length(seed_gene_indexs), ncol = length(genes) / nrow(genes));
  non_seed_genes <- genes;
  for (seed_gene_index_i in seq_len(length(seed_gene_indexs))) {
    seed_genes[seed_gene_index_i, ] <-
      genes[(seed_gene_indexs[seed_gene_index_i]), ];
  };
  non_seed_genes <- genes[-seed_gene_indexs, ];
  return(list(seed_genes, non_seed_genes));
};
