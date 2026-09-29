# =============================================================================
# 06  AIM 3 (descriptive) - does the sex difference change across a menopause proxy?
# Donors < 50 y vs >= 50 y (primary; only 2 women were < 45 in the inventory); sensitivity: < 45 vs >= 55.
# Model per dataset: y (standardized) ~ sex * older [+ assay]; interaction = change in (women - men)
# from the younger to the older band. A dataset contributes only with >= 2 women and >= 2 men per band.
# Reported with CIs only; this aim is not powered for inference.
# Makes: aim3.rds
# =============================================================================
A1 <- readRDS("aim1.rds"); A2 <- readRDS("aim2.rds")
dat <- A1$scores |> left_join(A2$S |> select(sample, peri_cap_ratio), by = "sample")
vars <- c(PRIMARY, "peri_cap_ratio")

band_fit <- function(d, y, cut_lo, cut_hi) {
  # donors known only by decade ("sixth decade" = 50-59, coded 55) straddle the excluded 45-54 band
  d <- d |> filter(is.finite(.data[[y]]), age < cut_lo | age >= cut_hi,
                  !(age_decade & cut_lo == 45 & age == 55)) |> mutate(older = age >= cut_hi)
  ok <- d |> count(sex, older)
  if (nrow(ok) < 4 || any(ok$n < 2)) return(NULL)
  d <- droplevels(d); d$yy <- as.numeric(scale(d[[y]]))
  f <- if (n_distinct(d$assay) > 1) yy ~ sex * older + assay else yy ~ sex * older
  co <- summary(lm(f, data = d))$coefficients
  term <- "sexfemale:olderTRUE"
  if (!term %in% rownames(co)) return(NULL)
  tibble(est = co[term, 1], se = co[term, 2], p = co[term, 4],
         n_w = sum(d$sex == "female"), n_m = sum(d$sex == "male"),
         young_w = sum(d$sex == "female" & !d$older), old_w = sum(d$sex == "female" & d$older))
}
run_band <- function(cut_lo, cut_hi) bind_rows(lapply(vars, function(v) {
  per <- bind_rows(lapply(split(dat, dat$dataset_id), function(d) {
    r <- band_fit(d, v, cut_lo, cut_hi); if (is.null(r)) NULL else mutate(r, dataset_id = d$dataset_id[1])
  }))
  if (!nrow(per)) return(NULL)
  meta_one(per) |> mutate(program = v)
}))

cat("\n== Donors by sex and age band ==\n")
dat |> mutate(band = cut(age, c(17, 44.99, 54.99, 120), labels = c("<45", "45-54", "55+"))) |>
  count(sex, band) |> pivot_wider(names_from = band, values_from = n, values_fill = 0) |> print()

main <- run_band(50, 50) |> mutate(cut = "<50 vs >=50")
sens <- run_band(45, 55) |> mutate(cut = "<45 vs >=55")
cat("\n== AIM 3: change in (women - men) from younger to older band, SD units ==\n")
bind_rows(main, sens) |> select(cut, program, k, n_w, n_m, est, lo, hi, p, I2) |>
  mutate(across(where(is.double), ~ signif(.x, 3))) |> print(n = Inf, width = Inf)

saveRDS(list(main = main, sens = sens), "aim3.rds")
