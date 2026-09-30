# =============================================================================
# 11  EXPLORATORY - which autosomal genes drive the donor clusters, and why?
# Run after 10_clustering_exploratory.R (uses lc, S, auto from that script).
# Splits the 2,000 most variable autosomal genes into 4 gene modules (same clustering as the heatmap rows),
# lists each module's top genes, checks them against contamination / technical signatures,
# and relates each donor's module score to sex, age, dataset, cell numbers and stress.
# Makes: tables/TableS_cluster_modules_genes.csv, tables/TableS_cluster_modules_donors.csv
# =============================================================================
m  <- lc[auto, ]; m <- m[head(order(apply(m, 1, var), decreasing = TRUE), 2000), ]
z  <- t(scale(t(m)))
gm <- cutree(hclust(as.dist(1 - cor(t(z))), method = "average"), k = 4)

sig <- list(
  cardiomyocyte_ambient = c("MYH7", "MYH6", "MYL2", "MYL3", "TNNT2", "TNNI3", "ACTC1", "MB", "TTN", "RYR2", "NPPA", "NPPB", "CKM", "COX6A2", "MYL7"),
  immune = c("PTPRC", "CD74", "HLA-DRA", "LYZ", "CD163", "C1QA", "C1QB", "F13A1", "CD14"),
  mural = c("RGS5", "ACTA2", "TAGLN", "MYH11", "PDGFRB", "NOTCH3", "KCNJ8"),
  fibroblast = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "PDGFRA"),
  stress = STRESS,
  mitochondrial = grep("^MT-", rownames(lc), value = TRUE),
  ribosomal = grep("^RP[SL]", rownames(lc), value = TRUE))

cat("\n== Gene modules: size, top genes (by variance), overlap with known signatures ==\n")
mods <- lapply(sort(unique(gm)), function(k) {
  g <- names(gm)[gm == k]
  top <- head(g[order(apply(m[g, , drop = FALSE], 1, var), decreasing = TRUE)], 25)
  ov <- sapply(sig, function(s) sum(g %in% s))
  cat("\nModule", k, "(", length(g), "genes)\n  top:", paste(top, collapse = ", "), "\n  overlap:",
      paste(names(ov)[ov > 0], ov[ov > 0], sep = "=", collapse = "; "), "\n")
  tibble(module = k, gene = g)
})
genes_tab <- bind_rows(mods)

# donor-level module scores and what they track
score <- sapply(sort(unique(gm)), function(k) colMeans(z[gm == k, , drop = FALSE]))
colnames(score) <- paste0("module", sort(unique(gm)))
stress <- colMeans(t(scale(t(lc[intersect(STRESS, rownames(lc)), ]))), na.rm = TRUE)
cm <- colMeans(t(scale(t(lc[intersect(sig$cardiomyocyte_ambient, rownames(lc)), ]))), na.rm = TRUE)
don <- S |> select(sample, donor_id, dataset, sex, age, n_capillary, n_cells) |>
  mutate(stress_score = stress, cardiomyocyte_ambient = cm) |> bind_cols(as_tibble(score))

cat("\n== Correlation (Spearman) of module scores with donor features ==\n")
feat <- c("age", "n_capillary", "n_cells", "stress_score", "cardiomyocyte_ambient")
print(round(sapply(colnames(score), function(k) sapply(feat, function(f)
  cor(don[[k]], don[[f]], method = "spearman", use = "complete.obs"))), 2))
cat("\n== Module scores: female minus male (Wilcoxon p) ==\n")
print(sapply(colnames(score), function(k) signif(wilcox.test(don[[k]] ~ don$sex)$p.value, 2)))

cat("\n== Donors with the most extreme module scores (|score| > 1) ==\n")
don |> filter(if_any(starts_with("module"), ~ abs(.x) > 1)) |>
  mutate(across(where(is.double), ~ round(.x, 2))) |> arrange(dataset) |> print(n = Inf, width = Inf)

write.csv(genes_tab, "tables/TableS_cluster_modules_genes.csv", row.names = FALSE)
write.csv(don, "tables/TableS_cluster_modules_donors.csv", row.names = FALSE)
