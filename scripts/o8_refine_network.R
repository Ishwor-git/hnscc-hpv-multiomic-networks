# 08_refine_networks.R — see docs/08_refine_networks.md

library(igraph)
library(dplyr)
library(Hmisc)

wgcna <- readRDS("data/processed/wgcna_results.rds")
datExpr <- wgcna$datExpr
traits <- wgcna$trait_data  # column is "hpv_status" (lowercase)

# datExpr's columns are Ensembl IDs (07 ran WGCNA directly on log_cpm), but
# deg_list/dmg_list use gene symbols - map here so DEG/DMG matching works.
# Keep the highest-MAD version of any symbol that maps from multiple
# Ensembl IDs (datExpr columns are already MAD-sorted, so !duplicated()
# keeps the first/highest-variance occurrence).
expr_raw <- readRDS("data/raw/tcga/expression/expression_raw.rds")
gene_map <- as.data.frame(SummarizedExperiment::rowData(expr_raw)) %>%
  select(gene_id, gene_name) %>% distinct()
rm(expr_raw); gc()

symbol_for_col <- gene_map$gene_name[match(colnames(datExpr), gene_map$gene_id)]
valid <- !is.na(symbol_for_col) & symbol_for_col != "" & !duplicated(symbol_for_col)

datExpr <- datExpr[, valid]
colnames(datExpr) <- symbol_for_col[valid]
module_colors <- setNames(wgcna$moduleColors[valid], colnames(datExpr))

cat("Genes retained after Ensembl->symbol mapping:", ncol(datExpr), "\n")

deg_list <- read.csv("results/tables/deg_list.csv")
dmg_list <- read.csv("results/tables/dmg_list.csv")

sig_module_names <- gsub("^ME", "", wgcna$sig_modules)
cat("Significant modules to refine:", paste(sig_module_names, collapse = ", "), "\n")

# optional TF list (dorothea high-confidence regulons as a TFcheckpoint substitute)
tf_genes <- tryCatch({
  if (!requireNamespace("dorothea", quietly = TRUE)) stop("dorothea not installed")
  unique(dorothea::dorothea_hs$tf[dorothea::dorothea_hs$confidence %in% c("A", "B")])
}, error = function(e) {
  cat("TF list unavailable (dorothea not installed) - TF annotation skipped.\n")
  character(0)
})

# optional PPI via STRINGdb
get_string_edges <- function(genes) {
  tryCatch({
    if (!requireNamespace("STRINGdb", quietly = TRUE)) stop("STRINGdb not installed")
    string_db <- STRINGdb::STRINGdb$new(version = "12.0", species = 9606,
                                        score_threshold = 700)
    mapped <- string_db$map(data.frame(gene = genes), "gene",
                            removeUnmappedRows = TRUE)
    ints <- string_db$get_interactions(mapped$STRING_id)
    ints %>%
      left_join(mapped, by = c("from" = "STRING_id")) %>%
      rename(gene1 = gene) %>%
      left_join(mapped, by = c("to" = "STRING_id")) %>%
      rename(gene2 = gene) %>%
      select(gene1, gene2) %>%
      filter(!is.na(gene1), !is.na(gene2))
  }, error = function(e) {
    cat("PPI lookup unavailable (STRINGdb not installed/reachable) - PPI edges skipped.\n")
    data.frame(gene1 = character(0), gene2 = character(0))
  })
}

# mutation data for significance annotation
mutation_data <- tryCatch(readRDS("data/processed/mutation_normalized.rds"),
                          error = function(e) NULL)

fisher_mutated_genes <- function(genes, hpv_status_vec) {
  if (is.null(mutation_data)) return(character(0))
  mb <- mutation_data$mutation_binary
  genes_present <- intersect(genes, colnames(mb))
  if (length(genes_present) == 0) return(character(0))
  sig <- c()
  for (g in genes_present) {
    tab <- table(mb[[g]][match(names(hpv_status_vec), mb$patient_barcode)],
                 hpv_status_vec)
    if (all(dim(tab) == c(2, 2))) {
      p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA)
      if (!is.na(p) && p <= 0.05) sig <- c(sig, g)
    }
  }
  sig
}

build_hpv_network <- function(module_name, hpv_label) {
  genes_in_module <- names(module_colors)[module_colors == module_name]
  sample_idx <- which(traits$hpv_status == ifelse(hpv_label == "HPV+", 1, 0))
  if (length(sample_idx) < 5 || length(genes_in_module) < 2) return(NULL)
  
  expr_subset <- datExpr[sample_idx, genes_in_module, drop = FALSE]
  cor_res <- rcorr(expr_subset, type = "spearman")
  
  r_mat <- cor_res$r
  p_mat <- cor_res$P
  edges <- which(abs(r_mat) >= 0.65 & p_mat <= 0.01, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]
  
  if (nrow(edges) == 0) return(NULL)
  
  edge_df <- data.frame(
    gene1 = rownames(r_mat)[edges[, 1]],
    gene2 = colnames(r_mat)[edges[, 2]],
    rho = r_mat[edges]
  )
  
  deg_genes_here <- intersect(genes_in_module, deg_list$gene_symbol)
  connected_to_deg <- unique(c(
    edge_df$gene1[edge_df$gene2 %in% deg_genes_here],
    edge_df$gene2[edge_df$gene1 %in% deg_genes_here]
  ))
  keep_genes <- union(deg_genes_here, connected_to_deg)
  
  edge_df <- edge_df %>% filter(gene1 %in% keep_genes, gene2 %in% keep_genes)
  if (nrow(edge_df) == 0) return(NULL)
  
  string_edges <- get_string_edges(keep_genes)
  if (nrow(string_edges) > 0) {
    string_edges <- string_edges %>%
      filter(gene1 %in% keep_genes, gene2 %in% keep_genes) %>%
      mutate(rho = NA_real_, edge_type = "PPI")
    edge_df$edge_type <- "co-expression"
    edge_df <- bind_rows(edge_df, string_edges)
  } else {
    edge_df$edge_type <- "co-expression"
  }
  
  g <- graph_from_data_frame(edge_df, directed = FALSE)
  
  hpv_vec <- setNames(traits$hpv_status, rownames(datExpr))
  mutated_sig <- fisher_mutated_genes(V(g)$name, hpv_vec)
  
  V(g)$is_DEG <- V(g)$name %in% deg_genes_here
  V(g)$is_DMG <- V(g)$name %in% dmg_list$gene_symbol
  V(g)$is_TF  <- V(g)$name %in% tf_genes
  V(g)$is_significantly_mutated <- V(g)$name %in% mutated_sig
  
  g
}

dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)
connection_metrics <- list()

for (mod in sig_module_names) {
  for (hpv in c("HPV+", "HPV-")) {
    cat("Building network:", mod, hpv, "\n")
    g <- build_hpv_network(mod, hpv)
    if (is.null(g)) {
      cat("  No network produced (insufficient edges/genes).\n")
      next
    }
    fname <- paste0("results/networks/", mod, "_",
                    ifelse(hpv == "HPV+", "hpv_pos", "hpv_neg"), ".graphml")
    write_graph(g, fname, format = "graphml")
    connection_metrics[[paste(mod, hpv)]] <- data.frame(
      module = mod, hpv_status = hpv,
      nodes = vcount(g), edges = ecount(g)
    )
    cat("  Saved:", fname, "| nodes:", vcount(g), "edges:", ecount(g), "\n")
  }
}

metrics_df <- bind_rows(connection_metrics)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(metrics_df, "results/tables/network_connection_metrics.csv", row.names = FALSE)

cat("\nSaved network connection metrics to results/tables/network_connection_metrics.csv\n")
print(metrics_df)