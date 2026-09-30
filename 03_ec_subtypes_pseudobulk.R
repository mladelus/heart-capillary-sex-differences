# =============================================================================
# 03  EC SUBTYPES + PSEUDOBULKS
# Atlases label endothelial subtypes differently, so every blood EC is re-assigned with the same rule:
# mean log(1 + CP10k) of capillary, arterial and venous markers, z-scored within dataset; the cell
# takes the subtype with the highest z. Agreement with the atlas labels is printed as a check.
# Then one pseudobulk per donor for capillary / arterial / venous ECs and for every other cell class.
# Makes: pb.rds (list: samples = donor table, mats = list of gene x donor matrices per cell group)
# =============================================================================
donors <- readRDS("inv_donors.rds") |> filter(included)
files  <- list.files("dl", full.names = TRUE)
dl     <- lapply(files, readRDS)
dl     <- Filter(function(x) any(donors$dataset_id == x$dataset_id & donors$donor_id == x$donor_id), dl)
message(length(dl), " donors loaded")

# ---- 1. Marker scores for every blood EC ----
ec_tab <- bind_rows(lapply(dl, function(x) {
  m <- x$ec_counts
  if (ncol(m) == 0) return(NULL)
  lib <- Matrix::colSums(m)
  sc <- sapply(ec_markers, function(g) {
    g <- intersect(g, rownames(m))
    if (!length(g)) return(rep(0, ncol(m)))
    Matrix::colMeans(log1p(t(t(m[g, , drop = FALSE]) / lib) * 1e4))
  })
  if (is.null(dim(sc))) sc <- matrix(sc, nrow = 1, dimnames = list(NULL, names(ec_markers)))
  tibble(dataset_id = x$dataset_id, donor_id = x$donor_id, soma_joinid = colnames(m),
         cell_type = x$meta$cell_type[match(colnames(m), as.character(x$meta$soma_joinid))],
         umis = lib) |> bind_cols(as_tibble(sc))
}))

ec_tab <- ec_tab |> group_by(dataset_id) |>
  mutate(across(c(capillary, arterial, venous), ~ { z <- as.numeric(scale(.x)); z[!is.finite(z)] <- 0; z },
                .names = "z_{.col}")) |>
  ungroup() |>
  mutate(subtype = c("capillary", "arterial", "venous")[max.col(cbind(z_capillary, z_arterial, z_venous),
                                                                ties.method = "first")])

ds_names <- readRDS("inv_datasets.rds") |> select(dataset_id, dataset_title)
cat("\n== Assigned subtype vs atlas label (row %), by dataset ==\n")
ec_tab |> left_join(ds_names, by = "dataset_id") |>
  count(dataset_title, cell_type, subtype) |>
  group_by(dataset_title, cell_type) |> mutate(pct = round(100 * n / sum(n))) |> ungroup() |>
  select(-n) |> pivot_wider(names_from = subtype, values_from = pct, values_fill = 0) |>
  print(n = Inf, width = Inf)

# Marker check: capillary-assigned cells should be high in CA4/RGCC and low in GJA5/ACKR1
cat("\n== Mean marker score by assigned subtype ==\n")
ec_tab |> group_by(subtype) |>
  summarise(cells = n(), capillary = mean(capillary), arterial = mean(arterial), venous = mean(venous)) |>
  print()

# ---- 2. Pseudobulks ----
stack_vecs <- function(vl) {                    # named count vectors -> genes x samples matrix (0-filled)
  genes <- sort(unique(unlist(lapply(vl, names))))
  out <- sapply(vl, function(v) { z <- setNames(numeric(length(genes)), genes); z[names(v)] <- v; z })
  if (is.null(dim(out))) out <- matrix(out, ncol = length(vl), dimnames = list(genes, names(vl)))
  out
}

sid <- function(x) paste(substr(x$dataset_id, 1, 8), x$donor_id, sep = "|")
ec_pb <- list(capillary = list(), arterial = list(), venous = list())
other_pb <- list()
n_tab <- list()
for (x in dl) {
  s <- sid(x)
  a <- ec_tab |> filter(dataset_id == x$dataset_id, donor_id == x$donor_id)
  for (st in names(ec_pb)) {
    cols <- a$soma_joinid[which(a$subtype == st)]
    if (length(cols) >= MIN_CT) ec_pb[[st]][[s]] <- Matrix::rowSums(x$ec_counts[, cols, drop = FALSE])
  }
  for (cl in setdiff(names(x$pb), c("blood EC", "other")))
    if (x$n_class[[cl]] >= MIN_CT) other_pb[[cl]][[s]] <- x$pb[[cl]]
  n_tab[[s]] <- tibble(sample = s, dataset_id = x$dataset_id, donor_id = x$donor_id,
                       n_cells = sum(x$n_class),
                       n_EC = sum(nrow(a)), n_capillary = sum(a$subtype == "capillary"),
                       n_arterial = sum(a$subtype == "arterial"), n_venous = sum(a$subtype == "venous"),
                       n_pericyte = if ("pericyte" %in% names(x$n_class)) x$n_class[["pericyte"]] else 0L,
                       n_cardiomyocyte = if ("cardiomyocyte" %in% names(x$n_class)) x$n_class[["cardiomyocyte"]] else 0L,
                       n_fibroblast = if ("fibroblast" %in% names(x$n_class)) x$n_class[["fibroblast"]] else 0L)
}
mats <- c(lapply(ec_pb, stack_vecs), lapply(other_pb, stack_vecs))

samples <- bind_rows(n_tab) |>
  left_join(donors |> select(dataset_id, donor_id, dataset_title, sex, age, age_decade, assay, suspension), by = c("dataset_id", "donor_id")) |>
  mutate(sex = relevel(factor(sex), ref = "male"), assay = factor(assay))

cat("\n== Donors with >= ", MIN_CAP, " capillary ECs, by dataset and sex ==\n", sep = "")
samples |> filter(n_capillary >= MIN_CAP) |> count(dataset_title, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print(n = Inf, width = Inf)
cat("\n== Pseudobulks per cell group ==\n")
print(sapply(mats, ncol))

saveRDS(list(samples = samples, mats = mats), "pb.rds")
saveRDS(ec_tab, "ec_subtypes.rds")
rm(dl); gc(verbose = FALSE)
