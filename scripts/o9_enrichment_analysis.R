# 09_enrichment_analysis.R — see docs/09_enrichment_analysis.md

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(igraph)

network_files <- list.files("results/networks", pattern = "\\.graphml$", full.names = TRUE)
cat("Found", length(network_files), "network files to enrich.\n")

msigdb_hallmark <- tryCatch({
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr not installed")
  msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    select(gs_name, gene_symbol)
}, error = function(e) {
  cat("MSigDB hallmark sets unavailable (msigdbr not installed) - skipping H category.\n")
  NULL
})

run_enrichment_for_network <- function(filepath) {
  g <- read_graph(filepath, format = "graphml")
  genes <- V(g)$name
  label <- gsub("\\.graphml$", "", basename(filepath))

  cat("\nEnriching:", label, "(", length(genes), "genes )\n")

  entrez <- tryCatch(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(entrez) || nrow(entrez) == 0) {
    cat("  No genes mapped to Entrez IDs - skipping.\n")
    return(NULL)
  }

  go_res <- tryCatch({
    enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
              pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff = 1)
  }, error = function(e) NULL)

  kegg_res <- tryCatch({
    enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
               pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff = 1)
  }, error = function(e) NULL)

  hallmark_res <- if (!is.null(msigdb_hallmark)) {
    tryCatch(
      enricher(genes, TERM2GENE = msigdb_hallmark, pAdjustMethod = "fdr",
               pvalueCutoff = 1, qvalueCutoff = 1),
      error = function(e) NULL
    )
  } else NULL

  results_list <- list()
  if (!is.null(go_res))       results_list$GO_BP  <- as.data.frame(go_res) %>% mutate(category = "GO_BP")
  if (!is.null(kegg_res))     results_list$KEGG   <- as.data.frame(kegg_res) %>% mutate(category = "KEGG")
  if (!is.null(hallmark_res)) results_list$MSigDB <- as.data.frame(hallmark_res) %>% mutate(category = "MSigDB_H")

  if (length(results_list) == 0) return(NULL)

  combined <- bind_rows(results_list) %>%
    filter(p.adjust <= 0.05) %>%
    mutate(network = label)

  cat("  Significant terms (FDR<=0.05):", nrow(combined), "\n")
  combined
}

all_enrichment <- lapply(network_files, run_enrichment_for_network) %>% bind_rows()

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(all_enrichment, "results/tables/enrichment_results.csv", row.names = FALSE)

cat("\nTotal significant enrichment terms across all networks:", nrow(all_enrichment), "\n")
cat("Saved to results/tables/enrichment_results.csv\n")
