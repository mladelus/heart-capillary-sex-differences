# =============================================================================
# 01  INVENTORY - which healthy adult heart ventricles are in the Census? (metadata only)
# Selects eligible donors and datasets by the pre-specified rules in 00_setup.R.
# Makes: inv_cells.rds (cell metadata for eligible donors), inv_donors.rds, inv_datasets.rds
# =============================================================================
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat()) |>
  select(dataset_id, dataset_title, collection_name, collection_doi)

obs <- census$get("census_data")$get("homo_sapiens")$obs$read(
  value_filter = "tissue_general == 'heart' & disease == 'normal' & is_primary_data == TRUE",
  column_names = c("soma_joinid", "dataset_id", "donor_id", "sex", "development_stage", "tissue",
                   "cell_type", "assay", "suspension_type", "raw_sum"))$concat() |>
  as.data.frame() |> as_tibble() |>
  mutate(across(where(is.factor), as.character),
         age = parse_age(development_stage), age_decade = development_stage %in% names(dec_mid),
         class = cell_class(cell_type)) |>
  left_join(ds, by = "dataset_id")

# ---- 1. What heart regions exist (check the VENT pattern) ----
cat("\n== Heart tissue labels, adults (cells) ==\n")
obs |> filter(age >= 18) |> count(tissue, sort = TRUE) |>
  mutate(in_VENT = str_detect(tissue, VENT)) |> print(n = Inf)

# ---- 2. Eligible donors: adult, ventricular tissue, >= MIN_EC endothelial cells, known sex ----
vent <- obs |> filter(age >= 18, str_detect(tissue, VENT), sex %in% c("female", "male"))
donors <- vent |> group_by(dataset_id, dataset_title, collection_name, donor_id, sex, age, age_decade) |>
  summarise(cells = n(), ECs = sum(class == "blood EC"),
            census_cap = sum(cell_type == "capillary endothelial cell"),
            pericytes = sum(class == "pericyte"),
            assay = paste(sort(unique(assay)), collapse = "; "),
            suspension = paste(sort(unique(suspension_type)), collapse = "; "),
            regions = paste(sort(unique(tissue)), collapse = "; "), .groups = "drop") |>
  mutate(eligible = ECs >= MIN_EC)

datasets <- donors |> filter(eligible) |>
  group_by(dataset_id, dataset_title, collection_name) |>
  summarise(women = sum(sex == "female"), men = sum(sex == "male"),
            age_w = paste(range(age[sex == "female"]), collapse = "-"),
            age_m = paste(range(age[sex == "male"]), collapse = "-"),
            median_EC = median(ECs), cells = sum(cells),
            assays = paste(sort(unique(assay)), collapse = " | "),
            suspension = paste(sort(unique(suspension)), collapse = " | "), .groups = "drop") |>
  mutate(included = women >= MIN_PER_SEX & men >= MIN_PER_SEX) |>
  left_join(ds |> select(dataset_id, collection_doi), by = "dataset_id") |>
  arrange(desc(women + men))

cat("\n== Datasets with eligible donors (>= ", MIN_EC, " ventricular ECs) ==\n", sep = "")
datasets |> select(-dataset_id) |> print(n = Inf, width = Inf)

inc <- datasets |> filter(included)
donors <- donors |> mutate(included = eligible & dataset_id %in% inc$dataset_id)

# ---- 3. Sex x assay balance within included datasets (confounding check) ----
cat("\n== Included donors by dataset x assay x sex ==\n")
donors |> filter(included) |> count(dataset_title, assay, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print(n = Inf, width = Inf)

# ---- 4. Possible donor overlap across datasets (same collection, same sex and age) ----
cat("\n== Donor IDs seen in more than one included dataset ==\n")
donors |> filter(included) |> group_by(donor_id) |> filter(n_distinct(dataset_id) > 1) |>
  select(donor_id, sex, age, dataset_title) |> print(n = Inf)
cat("\n== Included datasets from the same collection (check for shared donors) ==\n")
inc |> count(collection_name) |> filter(n > 1) |> print()

# ---- 5. Age bands (menopause proxy) ----
cat("\n== Included donors by sex and age band ==\n")
donors |> filter(included) |>
  mutate(band = cut(age, c(17, 44, 54, 120), labels = c("18-44", "45-54", "55+"))) |>
  count(sex, band) |> pivot_wider(names_from = band, values_from = n, values_fill = 0) |> print()

cat("\nTotal included: ", sum(donors$included), " donors (",
    sum(donors$included & donors$sex == "female"), " women, ",
    sum(donors$included & donors$sex == "male"), " men) in ", nrow(inc), " datasets; ",
    sum(vent$donor_id %in% donors$donor_id[donors$included] &
        vent$dataset_id %in% inc$dataset_id), " cells to download.\n", sep = "")

saveRDS(vent |> semi_join(donors |> filter(included), by = c("dataset_id", "donor_id")), "inv_cells.rds")
saveRDS(donors, "inv_donors.rds")
saveRDS(datasets, "inv_datasets.rds")
