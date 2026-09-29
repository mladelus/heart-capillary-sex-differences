# =============================================================================
# 07  AIM 4 - is the X-Y paralog pattern in capillary ECs specific to capillaries?
# Same donors, several cell groups (capillary ECs, arterial ECs, cardiomyocytes, fibroblasts,
# pericytes, myeloid cells, smooth muscle). For each X-Y pair:
#   X copy alone and combined X + Y (log2 CPM), women - men, per cell group (lm ~ sex + age [+ assay])
#   Y share in men = Y CPM / (X + Y CPM)
#   specificity: lmer(expr ~ sex * group + age [+ assay] + (1 | donor)), capillary as reference;
#   sex:group terms = (women - men in that group) - (women - men in capillaries); meta across datasets.
# X and Y copies are also reported separately, because the paralogs are not functionally equivalent.
# Makes: aim4.rds
# =============================================================================
P <- readRDS("pb.rds")
groups <- intersect(c("capillary", "arterial", "cardiomyocyte", "fibroblast", "pericyte", "myeloid", "smooth muscle"),
                    names(P$mats))
xy_genes <- c(xy_pairs$X, xy_pairs$Y, "XIST", "JPX")

long <- bind_rows(lapply(groups, function(g) {
  m <- P$mats[[g]]
  S <- P$samples |> filter(sample %in% colnames(m))
  bind_rows(lapply(split(S, S$dataset_id), function(d) {
    cp <- edgeR::cpm(calcNormFactors(DGEList(m[, d$sample, drop = FALSE])))      # CPM, all genes as library
    gg <- intersect(xy_genes, rownames(cp))
    as_tibble(t(cp[gg, , drop = FALSE])) |> mutate(sample = d$sample) |>
      pivot_longer(-sample, names_to = "gene", values_to = "cpm") |>
      left_join(d |> select(sample, dataset_id, donor_id, sex, age, assay), by = "sample") |>
      mutate(group = g)
  }))
}))

# X alone, Y alone, X + Y per pair
pair_long <- bind_rows(lapply(seq_len(nrow(xy_pairs)), function(i) {
  x <- long |> filter(gene == xy_pairs$X[i]) |> select(sample, group, dataset_id, donor_id, sex, age, assay, X = cpm)
  y <- long |> filter(gene == xy_pairs$Y[i]) |> select(sample, group, Y = cpm)
  inner_join(x, y, by = c("sample", "group")) |>
    mutate(pair = paste(xy_pairs$X[i], xy_pairs$Y[i], sep = "/"),
           X_only = log2(X + 1), Y_only = log2(Y + 1), X_plus_Y = log2(X + Y + 1), Y_share = Y / (X + Y))
}))

# ---- 1. Women - men per group (log2 units), meta across datasets ----
lfit <- function(d, y) {
  d <- droplevels(d)
  if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX) return(NULL)
  f <- as.formula(paste(y, if (n_distinct(d$assay) > 1) "~ sex + age + assay" else "~ sex + age"))
  co <- summary(lm(f, data = d))$coefficients
  if (!"sexfemale" %in% rownames(co)) return(NULL)
  tibble(est = co["sexfemale", 1], se = co["sexfemale", 2], p = co["sexfemale", 4],
         n_w = sum(d$sex == "female"), n_m = sum(d$sex == "male"))
}
by_group <- pair_long |> pivot_longer(c(X_only, X_plus_Y), names_to = "measure", values_to = "value") |>
  group_by(pair, group, measure, dataset_id) |> group_modify(~ { r <- lfit(.x, "value"); if (is.null(r)) tibble() else r }) |>
  ungroup() |> group_by(pair, group, measure) |> group_modify(~ meta_one(.x)) |> ungroup() |>
  group_by(measure) |> mutate(FDR = p.adjust(p, "BH")) |> ungroup()

cat("\n== Women - men (log2), pooled across datasets; * FDR < 0.05 ==\n")
by_group |> mutate(v = paste0(sprintf("%+.2f", est), ifelse(FDR < 0.05, "*", ""))) |>
  select(measure, pair, group, v) |> pivot_wider(names_from = group, values_from = v) |>
  arrange(measure, pair) |> print(n = Inf, width = Inf)

cat("\n== XIST / JPX check (women - men, log2) ==\n")
long |> filter(gene %in% c("XIST", "JPX")) |> mutate(value = log2(cpm + 1)) |>
  group_by(gene, group, dataset_id) |> group_modify(~ { r <- lfit(.x, "value"); if (is.null(r)) tibble() else r }) |>
  ungroup() |> group_by(gene, group) |> group_modify(~ meta_one(.x)) |> ungroup() |>
  select(gene, group, k, est, lo, hi) |> mutate(across(where(is.double), ~ round(.x, 2))) |> print(n = Inf)

# ---- 2. Y share in men ----
yshare <- pair_long |> filter(sex == "male") |> group_by(pair, group) |>
  summarise(Y_share_pct = round(100 * median(Y_share, na.rm = TRUE)), .groups = "drop")
cat("\n== Y copy share in men (median %, all datasets) ==\n")
yshare |> pivot_wider(names_from = group, values_from = Y_share_pct) |> print(width = Inf)

# ---- 3. Is the capillary sex difference different from other cell groups? ----
spec <- pair_long |> pivot_longer(c(X_only, Y_only, X_plus_Y), names_to = "measure", values_to = "value") |>
  mutate(group = relevel(factor(group), ref = "capillary")) |>
  group_by(pair, measure, dataset_id) |>
  group_modify(~ {
    d <- droplevels(.x)
    if (!"capillary" %in% d$group || n_distinct(d$group) < 2 ||
        sum(d$sex[d$group == "capillary"] == "female") < MIN_PER_SEX ||
        sum(d$sex[d$group == "capillary"] == "male") < MIN_PER_SEX) return(tibble())
    f <- if (n_distinct(d$assay) > 1) value ~ sex * group + age + assay + (1 | donor_id) else value ~ sex * group + age + (1 | donor_id)
    co <- tryCatch(summary(suppressMessages(lmerTest::lmer(f, data = d)))$coefficients, error = function(e) NULL)
    if (is.null(co)) return(tibble())
    k <- grep("^sexfemale:group", rownames(co))
    tibble(group = sub("^sexfemale:group", "", rownames(co)[k]), est = co[k, 1], se = co[k, 2], p = co[k, ncol(co)],
           n_w = sum(d$sex == "female" & d$group == "capillary"), n_m = sum(d$sex == "male" & d$group == "capillary"))
  }) |> ungroup() |>
  group_by(pair, measure, group) |> group_modify(~ meta_one(.x)) |> ungroup() |>
  group_by(measure) |> mutate(FDR = p.adjust(p, "BH")) |> ungroup()

cat("\n== Specificity: (women - men in group) minus (women - men in capillary ECs), log2; * FDR < 0.05 ==\n")
spec |> mutate(v = paste0(sprintf("%+.2f", est), ifelse(FDR < 0.05, "*", ""))) |>
  select(measure, pair, group, v) |> pivot_wider(names_from = group, values_from = v) |>
  arrange(measure, pair) |> print(n = Inf, width = Inf)

saveRDS(list(long = long, pair_long = pair_long, by_group = by_group, yshare = yshare, spec = spec), "aim4.rds")
