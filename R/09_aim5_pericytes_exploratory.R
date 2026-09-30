# =============================================================================
# 09  AIM 5 (EXPLORATORY, added after the primary results) - pericytes
# Added because women had fewer pericytes per capillary EC (Aim 2, p = 0.12). Not pre-specified:
# report as exploratory, with estimates and CIs; Holm within the family is shown for reference only.
#  (a) pericyte function programs, women - men (same scoring, model and meta-analysis as Aim 1)
#  (b) pericyte-to-capillary ratio: per dataset, + log total cells, + age, and sex x age (< 50 vs >= 50)
# Makes: aim5.rds, tables/TableS_aim5_*.csv, figures/Figure6_pericytes.(pdf|tiff)
# =============================================================================
set.seed(2027)
P  <- readRDS("pb.rds")
pe <- P$mats$pericyte
S  <- P$samples |> filter(n_pericyte >= MIN_CT, sample %in% colnames(pe))

peri_programs <- list(
  Pericyte_identity = c("RGS5", "NOTCH3", "PDGFRB", "CSPG4", "HIGD1B"),
  Contractile       = c("ACTA2", "TAGLN", "MYL9", "CNN1", "MYH11"),
  KATP_channel      = c("KCNJ8", "ABCC9"),
  Constrictor_receptors = c("EDNRA", "AGTR1"),
  NO_cGMP_response  = c("GUCY1A1", "GUCY1B1", "PRKG1", "PDE5A"))
PERI <- names(peri_programs)
sets <- c(peri_programs, pos_control)

# ---- (a) programs ----
by_ds <- split(S, S$dataset_id)
by_ds <- by_ds[sapply(by_ds, function(d) sum(d$sex == "female") >= MIN_PER_SEX && sum(d$sex == "male") >= MIN_PER_SEX)]
lc_ds <- lapply(by_ds, function(d) {
  y <- DGEList(pe[, d$sample, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE])
  cpm(y, log = TRUE, prior.count = 1)
})
scores <- bind_rows(lapply(names(by_ds), function(k)
  bind_cols(by_ds[[k]], as_tibble(sapply(sets, function(g) prog_score(lc_ds[[k]], g))))))

cat("\n== Pericyte donors per dataset ==\n")
scores |> count(dataset_title, sex) |> pivot_wider(names_from = sex, values_from = n) |> print()
cat("\n== Program genes passing expression filter ==\n")
print(sapply(lc_ds, function(lc) sapply(sets, function(g) paste0(sum(g %in% rownames(lc)), "/", length(g)))))

fit_sets <- function(sc, progs, extra = NULL) bind_rows(lapply(progs, function(prog) bind_rows(lapply(
  split(sc, sc$dataset_id), function(d) {
    r <- sex_fit(d, prog, extra); if (is.null(r)) NULL else mutate(r, dataset_id = d$dataset_id[1], program = prog)
  }))))
meta_all <- function(per) per |> group_by(program) |> group_modify(~ meta_one(.x)) |> ungroup()

per_ds <- fit_sets(scores, names(sets))
res <- meta_all(per_ds) |>
  mutate(p_holm_ref = ifelse(program %in% PERI, p.adjust(ifelse(program %in% PERI, p, NA), "holm"), NA))
r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))
cat("\n== AIM 5a (exploratory): pericyte programs, women - men (SD units) ==\n")
res |> select(program, k, n_w, n_m, est, lo, hi, p, p_holm_ref, lo90, hi90, within_0.5, I2, same_dir) |>
  r3() |> print(width = Inf)
cat("\n== Per-dataset estimates ==\n")
per_ds |> left_join(S |> distinct(dataset_id, dataset_title), by = "dataset_id") |>
  select(program, dataset_title, n_w, n_m, est, se, p) |> r3() |> arrange(program) |> print(n = Inf, width = Inf)

# calibration: shuffle sex within dataset x assay
perm <- bind_rows(lapply(1:1000, function(i) {
  sc <- scores |> group_by(dataset_id, assay) |> mutate(sex = sample(sex)) |> ungroup()
  meta_all(fit_sets(sc, PERI)) |> select(program, p)
}))
calib <- perm |> group_by(program) |> summarise(false_pos_rate_at_0.05 = mean(p < 0.05), pp = list(p)) |>
  left_join(res |> select(program, p_obs = p), by = "program") |>
  rowwise() |> mutate(p_perm = (1 + sum(unlist(pp) <= p_obs)) / (1 + length(unlist(pp)))) |> ungroup() |> select(-pp)
cat("\n== Calibration (1,000 permutations) ==\n")
calib |> r3() |> print(width = Inf)

# ---- (b) pericyte-to-capillary ratio, robustness ----
R <- P$samples |> filter(n_EC >= MIN_EC) |>
  mutate(peri_cap_ratio = log((n_pericyte + 0.5) / (n_capillary + 0.5)), log_cells = log(n_cells),
         older = age >= 50)
ratio_fit <- function(d, f, term) {
  d <- droplevels(d)
  if (n_distinct(d$assay) > 1) f <- update(f, . ~ . + assay)
  co <- summary(lm(f, data = d |> mutate(yy = as.numeric(scale(peri_cap_ratio)))))$coefficients
  if (!term %in% rownames(co)) return(NULL)
  tibble(est = co[term, 1], se = co[term, 2], p = co[term, 4], n_w = sum(d$sex == "female"), n_m = sum(d$sex == "male"))
}
models <- list(
  "sex + age (Aim 2 model)" = list(yy ~ sex + age, "sexfemale"),
  "+ log total cells"       = list(yy ~ sex + age + log_cells, "sexfemale"),
  "sex x age band (< 50 vs >= 50)" = list(yy ~ sex * older, "sexfemale:olderTRUE"))
ratio <- bind_rows(lapply(names(models), function(m) {
  per <- bind_rows(lapply(split(R, R$dataset_id), function(d) {
    ok <- if (grepl("band", m)) all(table(d$sex, d$older) >= 2) && length(unique(d$older)) == 2 else TRUE
    if (!ok) return(NULL)
    r <- ratio_fit(d, models[[m]][[1]], models[[m]][[2]]); if (is.null(r)) NULL else mutate(r, dataset_id = d$dataset_id[1])
  }))
  bind_rows(per |> mutate(level = "dataset"), meta_one(per) |> mutate(level = "pooled")) |> mutate(model = m)
}))
cat("\n== AIM 5b: pericytes per capillary EC, women - men (SD units) ==\n")
ratio |> left_join(S |> distinct(dataset_id, dataset_title), by = "dataset_id") |>
  select(model, level, dataset_title, n_w, n_m, est, lo, hi, se, p) |> r3() |> print(n = Inf, width = Inf)
cat("\n== Raw pericytes per 100 capillary ECs, by dataset and sex (median, IQR) ==\n")
R |> group_by(dataset_title, sex) |>
  summarise(donors = n(), median = median(100 * n_pericyte / n_capillary),
            q1 = quantile(100 * n_pericyte / n_capillary, 0.25), q3 = quantile(100 * n_pericyte / n_capillary, 0.75),
            .groups = "drop") |> r3() |> print(width = Inf)

# ---- outputs ----
dir.create("tables", showWarnings = FALSE)
write.csv(res, "tables/TableS_aim5_pericyte_programs.csv", row.names = FALSE)
write.csv(per_ds, "tables/TableS_aim5_pericyte_programs_per_dataset.csv", row.names = FALSE)
write.csv(calib, "tables/TableS_aim5_calibration.csv", row.names = FALSE)
write.csv(ratio, "tables/TableS_aim5_pericyte_ratio.csv", row.names = FALSE)

lab5 <- c(Pericyte_identity = "Pericyte identity", Contractile = "Contractile", KATP_channel = "K-ATP channel",
          Constrictor_receptors = "Constrictor receptors", NO_cGMP_response = "NO / cGMP response",
          X_escape = "X-escape (positive control)")
m5 <- res |> mutate(label = factor(lab5[program], levels = rev(lab5)))
f6 <- ggplot(m5, aes(est, label)) +
  annotate("rect", xmin = -MARGIN, xmax = MARGIN, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.6) +
  geom_point(shape = 18, size = 3) +
  labs(x = "Female minus male in pericytes (donor-level SD units)", y = NULL) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave("figures/Figure6_pericytes.pdf", f6, width = 170, height = 75, units = "mm")
ggsave("figures/Figure6_pericytes.tiff", f6, width = 170, height = 75, units = "mm", dpi = 300, compression = "lzw")
saveRDS(list(scores = scores, per_ds = per_ds, res = res, calib = calib, ratio = ratio), "aim5.rds")
message("Aim 5 done.")
