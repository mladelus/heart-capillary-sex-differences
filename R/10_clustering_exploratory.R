# =============================================================================
# 10  EXPLORATORY - do donors cluster by sex? (unsupervised, capillary ECs)
# Capillary pseudobulks (45 donors) -> TMM log-CPM -> remove dataset/chemistry effects (limma::removeBatchEffect)
# -> 2,000 most variable genes -> hierarchical clustering (1 - Pearson correlation, average linkage).
# Two versions:
#   A. autosomal genes only (the real question: is there a genome-wide sex signature?)
#   B. all genes including X and Y (positive control: should separate by sex, driven by XIST/Y genes)
# Agreement between clusters and sex: adjusted Rand index (2-cluster cut) and nearest-neighbour
# same-sex rate, each against 1,000 sex-label permutations.
# Makes: figures/FigureS_clustering_autosomal.(pdf|png), figures/FigureS_clustering_all_genes.(pdf|png),
#        tables/TableS_clustering.csv
# =============================================================================
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
suppressPackageStartupMessages({ library(pheatmap); library(EnsDb.Hsapiens.v86) })
set.seed(2028)

P   <- readRDS("pb.rds")
cap <- P$mats$capillary
S   <- P$samples |> filter(n_capillary >= MIN_CAP, sample %in% colnames(cap)) |> droplevels()
S   <- S |> mutate(dataset = ifelse(grepl("^Heart$", dataset_title), "Harmonized cohort", "Heart Cell Atlas"))

y <- DGEList(cap[, S$sample])
y <- calcNormFactors(y[filterByExpr(y, group = S$sex), , keep.lib.sizes = FALSE])
lc <- cpm(y, log = TRUE, prior.count = 1)
# remove dataset and chemistry differences (sex is NOT given to the correction)
lc <- removeBatchEffect(lc, batch = S$dataset, batch2 = S$assay)

gg  <- genes(EnsDb.Hsapiens.v86, columns = c("gene_name", "seq_name"))
chr <- setNames(as.character(seqnames(gg)), gg$gene_name)
auto <- rownames(lc)[chr[rownames(lc)] %in% c(1:22)]

ann <- data.frame(Sex = as.character(S$sex), Dataset = S$dataset, Age = S$age, row.names = S$sample)
ann_col <- list(Sex = c(female = "#B2182B", male = "#2166AC"),
                Dataset = c("Harmonized cohort" = "grey30", "Heart Cell Atlas" = "grey75"))

cluster_sex <- function(genes, label, n_top = 2000) {
  m <- lc[genes, , drop = FALSE]
  top <- head(order(apply(m, 1, var), decreasing = TRUE), min(n_top, nrow(m)))
  m <- m[top, ]
  d <- as.dist(1 - cor(m))
  hc <- hclust(d, method = "average")
  sex <- as.character(S$sex)
  ari <- function(a, b) {                            # adjusted Rand index
    tab <- table(a, b); n <- sum(tab); s2 <- function(x) sum(choose(x, 2))
    e <- s2(rowSums(tab)) * s2(colSums(tab)) / choose(n, 2)
    (s2(tab) - e) / ((s2(rowSums(tab)) + s2(colSums(tab))) / 2 - e)
  }
  nn_same <- function(sx) { dm <- as.matrix(d); diag(dm) <- Inf; mean(sx == sx[apply(dm, 1, which.min)]) }
  cl2 <- cutree(hc, 2)
  obs <- c(ARI = ari(cl2, sex), NN_same_sex = nn_same(sex))
  perm <- replicate(1000, { s <- sample(sex); c(ari(cl2, s), nn_same(s)) })
  res <- tibble(genes = label, n_genes = nrow(m), ARI = obs[1], ARI_perm_p = mean(perm[1, ] >= obs[1]),
                NN_same_sex = obs[2], NN_perm_mean = mean(perm[2, ]), NN_perm_p = mean(perm[2, ] >= obs[2]))
  z <- t(scale(t(m)))
  z[z > 3] <- 3; z[z < -3] <- -3
  file <- file.path("figures", paste0("FigureS_clustering_", gsub(" ", "_", label)))
  for (ext in c("pdf", "png"))
    pheatmap(z, cluster_rows = TRUE, cluster_cols = hc, show_rownames = FALSE, show_colnames = FALSE,
             annotation_col = ann, annotation_colors = ann_col, treeheight_row = 0,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
             main = paste0("Capillary ECs, 45 donors: ", nrow(m), " most variable ", label, " genes"),
             filename = paste0(file, ".", ext), width = 7, height = 6)
  res
}

sexchr <- rownames(lc)[chr[rownames(lc)] %in% c("X", "Y")]
res <- bind_rows(cluster_sex(auto, "autosomal"), cluster_sex(rownames(lc), "all genes"),
                 cluster_sex(sexchr, "X and Y", n_top = 200))
cat("\n== Do capillary ECs cluster by sex? ==\n")
res |> mutate(across(where(is.double), ~ signif(.x, 3))) |> print(width = Inf)
cat("ARI: 0 = no agreement with sex, 1 = perfect. NN_same_sex: share of donors whose most similar donor is the same sex",
    "(expected by chance ~", round(mean(res$NN_perm_mean), 2), ").\n")
write.csv(res, "tables/TableS_clustering.csv", row.names = FALSE)

# ---- What does drive donor-to-donor variation? PCA of the 2,000 most variable autosomal genes ----
m  <- lc[auto, ]; m <- m[head(order(apply(m, 1, var), decreasing = TRUE), 2000), ]
pc <- prcomp(t(m), scale. = FALSE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
stress <- colMeans(t(scale(t(lc[intersect(STRESS, rownames(lc)), ]))), na.rm = TRUE)
cov <- data.frame(sex = S$sex, age = S$age, log_capillary_cells = log(S$n_capillary),
                  log_all_cells = log(S$n_cells), stress_score = stress)
r2 <- sapply(1:5, function(k) sapply(names(cov), function(v) summary(lm(pc$x[, k] ~ cov[[v]]))$r.squared))
colnames(r2) <- paste0("PC", 1:5, " (", ve[1:5], "%)")
cat("\n== Variance in each principal component explained by each factor (R^2) ==\n")
print(round(r2, 2))
write.csv(data.frame(factor = rownames(r2), round(r2, 3)), "tables/TableS_clustering_PC_drivers.csv", row.names = FALSE)
message("Heatmaps in ", file.path(getwd(), "figures"))
