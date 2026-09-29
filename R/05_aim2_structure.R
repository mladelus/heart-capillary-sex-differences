# =============================================================================
# 05  AIM 2 - microvascular structure proxies and capillary-pericyte signaling
# Per donor (ventricular cells):
#   pericyte-to-capillary-EC ratio  = log((pericytes + 0.5) / (capillary ECs + 0.5))
#   capillary share of blood ECs    = logit((capillary + 0.5) / (blood ECs + 1))
#   EC share of all cells           = logit((blood ECs + 0.5) / (all cells + 1))
# Composition depends on dissociation (cells vs nuclei), so results are also split by suspension type.
# Signaling: for each pair, donor score = z(ligand log-CPM in sender) + z(receptor log-CPM in receiver),
#   using capillary and pericyte pseudobulks from the same donor.
# Same model and meta-analysis as Aim 1 (women - men, SD units; Holm within each family).
# Makes: aim2.rds
# =============================================================================
P <- readRDS("pb.rds")
S <- P$samples |>
  mutate(peri_cap_ratio = log((n_pericyte + 0.5) / (n_capillary + 0.5)),
         cap_share      = qlogis((n_capillary + 0.5) / (n_EC + 1)),
         EC_share       = qlogis((n_EC + 0.5) / (n_cells + 1)))
comp_vars <- c("peri_cap_ratio", "cap_share", "EC_share")

fit_vars <- function(sc, vars) bind_rows(lapply(vars, function(v) bind_rows(lapply(split(sc, sc$dataset_id), function(d) {
  r <- sex_fit(d, v); if (is.null(r)) NULL else mutate(r, dataset_id = d$dataset_id[1], program = v)
}))))
meta_fam <- function(per) per |> group_by(program) |> group_modify(~ meta_one(.x)) |> ungroup() |>
  mutate(p_holm = p.adjust(p, "holm"), verdict = verdict(p_holm, lo90, hi90, same_dir, k))
r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))

# ---- 1. Composition ----
comp_per <- fit_vars(S |> filter(n_EC >= MIN_EC), comp_vars)
comp <- meta_fam(comp_per)
cat("\n== AIM 2a: composition, women - men (SD units) ==\n")
comp |> select(program, k, n_w, n_m, est, lo, hi, p, p_holm, lo90, hi90, I2, verdict) |> r3() |> print(width = Inf)
cat("\n== Composition by suspension type (cells vs nuclei) ==\n")
comp_susp <- bind_rows(lapply(split(comp_per |> left_join(S |> group_by(dataset_id) |> summarise(suspension = first(suspension)), by = "dataset_id"),
                                    ~ suspension), function(x) meta_fam(x |> select(-suspension)) |> mutate(suspension = x$suspension[1])))
comp_susp |> select(suspension, program, k, est, lo, hi, p) |> r3() |> print(width = Inf)

cat("\n== Raw medians by sex ==\n")
S |> filter(n_EC >= MIN_EC) |> group_by(sex) |>
  summarise(donors = n(), pericytes_per_100_cap = median(100 * n_pericyte / pmax(n_capillary, 1)),
            cap_pct_of_EC = median(100 * n_capillary / n_EC), EC_pct_of_cells = median(100 * n_EC / n_cells)) |>
  r3() |> print()

# ---- 2. Capillary-pericyte signaling ----
both <- intersect(colnames(P$mats$capillary), colnames(P$mats$pericyte))
Sb <- S |> filter(sample %in% both, n_capillary >= MIN_CAP, n_pericyte >= MIN_CT)
lc_of <- function(mat, cols) { y <- calcNormFactors(DGEList(mat[, cols, drop = FALSE])); cpm(y, log = TRUE, prior.count = 1) }
lr_scores <- bind_rows(lapply(Filter(function(d) sum(d$sex == "female") >= MIN_PER_SEX && sum(d$sex == "male") >= MIN_PER_SEX,
                                     split(Sb, Sb$dataset_id)), function(d) {
  lc <- list(capillary = lc_of(P$mats$capillary, d$sample), pericyte = lc_of(P$mats$pericyte, d$sample))
  zz <- function(v) if (is.null(v) || length(v) < 2 || !isTRUE(sd(v) > 0)) rep(NA_real_, nrow(d)) else as.numeric(scale(v))
  g  <- function(cell, gene) if (gene %in% rownames(lc[[cell]])) lc[[cell]][gene, ] else NULL
  sc <- sapply(seq_len(nrow(lr_pairs)), function(i) {
    a <- zz(g(lr_pairs$sender[i], lr_pairs$ligand[i])); b <- zz(g(lr_pairs$receiver[i], lr_pairs$receptor[i]))
    a + b
  })
  colnames(sc) <- lr_pairs$pair
  bind_cols(d, as_tibble(sc))
}))
lr_per <- fit_vars(lr_scores, lr_pairs$pair)
lr <- meta_fam(lr_per)
cat("\n== AIM 2b: capillary-pericyte signaling, women - men (SD units) ==\n")
cat("Donors with both pseudobulks: ", nrow(Sb), " (", sum(Sb$sex == "female"), " W, ", sum(Sb$sex == "male"), " M)\n", sep = "")
lr |> select(program, k, n_w, n_m, est, lo, hi, p, p_holm, lo90, hi90, I2, verdict) |> r3() |> print(width = Inf)

saveRDS(list(S = S, comp_per = comp_per, comp = comp, comp_susp = comp_susp,
             lr_scores = lr_scores, lr_per = lr_per, lr = lr), "aim2.rds")
