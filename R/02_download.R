# =============================================================================
# 02  DOWNLOAD - one donor at a time; keeps only what later steps need
# For each included donor (ventricular cells only):
#   - pseudobulk counts per broad cell class (all genes) and the number of cells per class
#   - the blood-EC count matrix (cells x genes, sparse) for uniform subtype assignment in 03
# Skips donors already saved, so if it stops (network, sleep), just run it again.
# Takes roughly 1-4 minutes per donor. Makes: dl/<dataset>__<donor>.rds
# =============================================================================
cells <- readRDS("inv_cells.rds")
dir.create("dl", showWarnings = FALSE)
keys <- cells |> distinct(dataset_id, donor_id)
message(nrow(keys), " donors to download")

for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  out <- file.path("dl", paste0(substr(k$dataset_id, 1, 8), "__", gsub("[^A-Za-z0-9-]", "_", k$donor_id), ".rds"))
  if (file.exists(out)) next
  message(format(Sys.time(), "%H:%M"), "  ", i, "/", nrow(keys), "  ", k$donor_id)
  m <- cells |> filter(dataset_id == k$dataset_id, donor_id == k$donor_id)
  seu <- get_seurat(census, organism = "Homo sapiens", obs_coords = m$soma_joinid,
                    obs_column_names = "soma_joinid")
  counts <- to_symbols(LayerData(seu, assay = "RNA", layer = "counts"))
  m <- m[match(seu$soma_joinid, m$soma_joinid), ]; rm(seu); gc(verbose = FALSE)

  n_class <- table(m$class)
  pb <- lapply(split(seq_len(ncol(counts)), m$class), function(j) Matrix::rowSums(counts[, j, drop = FALSE]))
  ec <- which(m$class == "blood EC")
  ec_counts <- counts[, ec, drop = FALSE]
  ec_counts <- ec_counts[Matrix::rowSums(ec_counts) > 0, , drop = FALSE]
  colnames(ec_counts) <- as.character(m$soma_joinid[ec])
  saveRDS(list(dataset_id = k$dataset_id, donor_id = k$donor_id,
               meta = m |> select(soma_joinid, cell_type, class, tissue, assay, raw_sum),
               n_class = n_class, pb = pb, ec_counts = ec_counts), out)
  rm(counts, ec_counts, pb); gc(verbose = FALSE)
}
message("Saved donors: ", length(list.files("dl")), " of ", nrow(keys))
