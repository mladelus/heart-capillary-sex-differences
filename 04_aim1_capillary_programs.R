# =============================================================================
# 04  AIM 1 (PRIMARY) - do women's and men's ventricular capillary ECs differ in function programs?
# Unit = donor. Per dataset: capillary pseudobulk -> TMM log-CPM -> program score (mean z of genes)
#   -> score standardized within dataset -> lm(score ~ sex + age [+ assay]) -> women - men in SD units.
# Datasets combined by random-effects meta-analysis (REML; fixed effect if < 3 datasets).
# Primary family: 6 programs, Holm-adjusted. Equivalence: 90% CI inside +/- 0.5 SD (TOST).
# Checks: X-escape positive control; stress-adjusted sensitivity; leave-one-dataset-out;
#         1,000 within-dataset sex-label permutations (false-positive calibration).
# Genome-wide: limma-voom per dataset, inverse-variance meta of gene-level sex effects.
# Makes: aim1.rds
# =============================================================================
set.seed(2026)
P <- readRDS("pb.rds")
cap <- P$mats$capillary
S <- P$samples |> filter(n_capillary >= MIN_CAP, sample %in% colnames(cap))

all_sets <- c(programs, pos_control, list(Stress = STRESS))

# ---- 1. Per-dataset scores ----
by_ds <- split(S, S$dataset_id)
by_ds <- by_ds[sapply(by_ds, function(d) sum(d$sex == "female") >= MIN_PER_SEX && sum(d$sex == "male") >= MIN_PER_SEX)]
cat("\n== Datasets entering Aim 1 ==\n")
bind_rows(lapply(by_ds, function(d) tibble(dataset = d$dataset_title[1], women = sum(d$sex == "female"),
                                           men = sum(d$sex == "male"), median_cap = median(d$n_capillary),
                                           assays = paste(unique(d$assay), collapse = " | ")))) |> print(width = Inf)

lc_ds <- lapply(by_ds, function(d) {
  y <- DGEList(cap[, d$sample, drop = FALSE])
  y <- y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  cpm(y, log = TRUE, prior.count = 1)
})

scores <- bind_rows(lapply(names(by_ds), function(k) {
  d <- by_ds[[k]]; lc <- lc_ds[[k]]
  sc <- sapply(all_sets, function(g) prog_score(lc, g))
  bind_cols(d, as_tibble(sc))
}))

cov_tab <- bind_rows(lapply(names(lc_ds), function(k) tibble(dataset_id = k,
  set = names(all_sets), genes_found = sapply(all_sets, function(g) sum(g %in% rownames(lc_ds[[k]]))),
  genes_total = lengths(all_sets))))
cat("\n== Program genes passing expression filter, per dataset ==\n")
cov_tab |> mutate(found = paste0(genes_found, "/", genes_total)) |> select(-genes_found, -genes_total) |>
  pivot_wider(names_from = dataset_id, values_from = found) |> print(width = Inf)

# ---- 2. Fits and meta-analysis ----
fit_all <- function(sc, sets, extra = NULL) {
  bind_rows(lapply(sets, function(p) bind_rows(lapply(split(sc, sc$dataset_id), function(d) {
    r <- sex_fit(d, p, extra); if (is.null(r)) NULL else mutate(r, dataset_id = d$dataset_id[1], program = p)
  }))))
}
run_meta <- function(per) per |> group_by(program) |> group_modify(~ meta_one(.x)) |> ungroup()

per_ds  <- fit_all(scores, names(all_sets))
main    <- run_meta(per_ds) |>
  mutate(p_holm = ifelse(program %in% PRIMARY, p.adjust(ifelse(program %in% PRIMARY, p, NA), "holm"), NA),
         MDE80 = (qnorm(0.975) + qnorm(0.8)) * (hi - lo) / (2 * qnorm(0.975)),
         verdict = ifelse(program %in% PRIMARY, verdict(p_holm, lo90, hi90, same_dir, k), NA))

xe <- main |> filter(program == "X_escape")
if (!nrow(xe) || !(xe$est > 0 && xe$p < 0.05)) warning("POSITIVE CONTROL FAILED: X-escape not higher in women")

r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))
cat("\n== AIM 1: women - men, SD units (random-effects meta) ==\n")
main |> select(program, k, n_w, n_m, est, lo, hi, p, p_holm, lo90, hi90, p_tost, within_0.5, pi_lo, pi_hi, I2,
               same_dir, MDE80, verdict) |> r3() |> print(n = Inf, width = Inf)
cat("\n== Per-dataset estimates ==\n")
per_ds |> left_join(readRDS("inv_datasets.rds") |> select(dataset_id, dataset_title), by = "dataset_id") |>
  select(program, dataset_title, n_w, n_m, est, se, p) |> r3() |> arrange(program) |> print(n = Inf, width = Inf)

# ---- 3. Sensitivity: + stress score, + log capillary cell number; leave-one-dataset-out ----
scores <- scores |> mutate(log_ncap = log(n_capillary))
sens <- bind_rows(
  run_meta(fit_all(scores, PRIMARY, "Stress")) |> mutate(analysis = "+ stress score"),
  run_meta(fit_all(scores, PRIMARY, "log_ncap")) |> mutate(analysis = "+ log capillary cells"))

# Sensitivity (added before outcome analysis, after the subtype check): capillaries as labeled by each atlas
cap_a <- P$mats$capillary_atlas
Sa <- P$samples |> filter(n_capillary_atlas >= MIN_CAP, sample %in% colnames(cap_a))
scores_a <- bind_rows(lapply(split(Sa, Sa$dataset_id), function(d) {
  if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX) return(NULL)
  y <- DGEList(cap_a[, d$sample, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE])
  lc <- cpm(y, log = TRUE, prior.count = 1)
  bind_cols(d, as_tibble(sapply(programs, function(g) prog_score(lc, g))))
}))
sens <- bind_rows(sens, run_meta(fit_all(scores_a, PRIMARY)) |> mutate(analysis = "atlas capillary labels"))
loo <- bind_rows(lapply(names(by_ds), function(k)
  run_meta(per_ds |> filter(dataset_id != k, program %in% PRIMARY)) |> mutate(left_out = k)))
cat("\n== Sensitivity analyses ==\n")
sens |> select(analysis, program, k, est, lo, hi, p) |> r3() |> print(n = Inf, width = Inf)
cat("\n== Leave one dataset out: range of pooled estimates ==\n")
loo |> group_by(program) |> summarise(min_est = min(est), max_est = max(est), min_p = min(p), max_p = max(p)) |>
  r3() |> print(width = Inf)

# ---- 4. Calibration: shuffle sex within dataset x assay, 1,000 times ----
perm_once <- function() {
  sc <- scores |> group_by(dataset_id, assay) |> mutate(sex = sample(sex)) |> ungroup()
  run_meta(fit_all(sc, c(PRIMARY, "X_escape"))) |> select(program, est, p)
}
perm <- bind_rows(lapply(1:1000, function(i) { if (i %% 100 == 0) message("permutation ", i); perm_once() }))
calib <- perm |> group_by(program) |>
  summarise(false_pos_rate_at_0.05 = mean(p < 0.05)) |>
  left_join(main |> select(program, p_obs = p), by = "program") |>
  left_join(perm |> group_by(program) |> summarise(pp = list(p)), by = "program") |>
  rowwise() |> mutate(p_perm = (1 + sum(unlist(pp) <= p_obs)) / (1 + length(unlist(pp)))) |> ungroup() |> select(-pp)
cat("\n== Calibration (1,000 sex-label permutations) ==\n")
calib |> r3() |> print(width = Inf)

# ---- 5. Genome-wide: limma-voom per dataset, inverse-variance meta ----
gw <- bind_rows(lapply(names(by_ds), function(k) {
  d <- droplevels(by_ds[[k]])
  y <- DGEList(cap[, d$sample, drop = FALSE]); y <- y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  X <- if (n_distinct(d$assay) > 1) model.matrix(~ sex + age + assay, d) else model.matrix(~ sex + age, d)
  fit <- eBayes(lmFit(voom(y, X), X))
  tt <- topTable(fit, coef = "sexfemale", number = Inf, sort.by = "none")
  tibble(dataset_id = k, gene = rownames(tt), logFC = tt$logFC, se = tt$logFC / tt$t, p = tt$P.Value)
}))
gmeta <- gw |> filter(is.finite(se), se > 0) |> group_by(gene) |> filter(n() >= 2) |>
  summarise(k = n(),
            same_dir = max(sum(logFC > 0), sum(logFC < 0)),
            pooled = sum(logFC / se^2) / sum(1 / se^2),
            Q_p = pchisq(sum(((logFC - pooled) / se)^2), n() - 1, lower.tail = FALSE),
            se = sqrt(1 / sum(1 / se^2)), .groups = "drop") |>
  rename(logFC = pooled) |>
  mutate(z = logFC / se, p = 2 * pnorm(-abs(z)), FDR = p.adjust(p, "BH"))
chrom <- tryCatch({
  suppressMessages(library(EnsDb.Hsapiens.v86))
  gg <- genes(EnsDb.Hsapiens.v86, columns = c("gene_name", "seq_name"))
  tibble(gene = gg$gene_name, chr = as.character(seqnames(gg))) |> distinct(gene, .keep_all = TRUE)
}, error = function(e) NULL)
if (!is.null(chrom)) gmeta <- gmeta |> left_join(chrom, by = "gene")

cat("\n== Genome-wide: positive-control genes ==\n")
gmeta |> filter(gene %in% c("XIST", "TSIX", "JPX", "KDM6A", "KDM5C", "RPS4Y1", "DDX3Y", "UTY", "KDM5D")) |>
  r3() |> print(width = Inf)
if ("chr" %in% names(gmeta)) {
  cat("\n== Genome-wide: FDR < 0.05 by chromosome class ==\n")
  gmeta |> mutate(chr_class = case_when(chr == "X" ~ "X", chr == "Y" ~ "Y", is.na(chr) ~ "unknown", TRUE ~ "autosome")) |>
    group_by(chr_class) |> summarise(tested = n(), FDR05 = sum(FDR < 0.05),
                                     FDR05_consistent = sum(FDR < 0.05 & same_dir == k & Q_p > 0.05)) |> print()
  cat("\n== Top autosomal genes (FDR < 0.05, same direction in every dataset) ==\n")
  gmeta |> filter(!chr %in% c("X", "Y"), FDR < 0.05, same_dir == k) |> arrange(p) |>
    select(gene, chr, k, logFC, se, p, FDR, Q_p) |> head(40) |> r3() |> print(n = 40, width = Inf)
}
# Replication of gene-level effects between datasets (Spearman of t-like z)
wide <- gw |> mutate(z = logFC / se) |> select(dataset_id, gene, z) |> pivot_wider(names_from = dataset_id, values_from = z)
if (ncol(wide) > 2) {
  cat("\n== Between-dataset correlation of gene-level sex effects (Spearman; autosomes only if known) ==\n")
  w <- if (!is.null(chrom)) wide |> left_join(chrom, by = "gene") |> filter(!chr %in% c("X", "Y")) |> select(-chr) else wide
  print(round(cor(as.matrix(w[, -1]), method = "spearman", use = "pairwise.complete.obs"), 3))
}

saveRDS(list(scores = scores, per_ds = per_ds, main = main, sens = sens, loo = loo, calib = calib,
             perm = perm, gw = gw, gmeta = gmeta, coverage = cov_tab), "aim1.rds")
