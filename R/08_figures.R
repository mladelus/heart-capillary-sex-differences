# =============================================================================
# 08  FIGURES (PDF + 300 dpi TIFF, 170 mm wide) and summary tables (CSV)
# Figure 2  Aim 1 forest plot: pooled women - men per program with equivalence bounds; dataset estimates
# Figure 3  Aim 1 checks: calibration (permutation p-values), leave-one-out, genome-wide volcano
# Figure 4  Aim 2 composition and capillary-pericyte signaling
# Figure 5  Aim 4 X-Y paralogs across cell groups
# =============================================================================
A1 <- readRDS("aim1.rds"); A2 <- readRDS("aim2.rds"); A3 <- readRDS("aim3.rds"); A4 <- readRDS("aim4.rds")
dsn <- readRDS("inv_datasets.rds") |> select(dataset_id, dataset_title)
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank()))
W <- "#B2182B"; M <- "#2166AC"
save_fig <- function(p, name, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = 170, height = h, units = "mm")
  ggsave(file.path("figures", paste0(name, ".tiff")), p, width = 170, height = h, units = "mm", dpi = 300, compression = "lzw")
}
lab <- c(NO_eNOS = "NO / eNOS", Endothelin_ACE = "Endothelin / ACE", Prostacyclin = "Prostacyclin",
         Barrier = "Barrier / junction", FA_transport = "Fatty-acid transport", Angiogenic_tip = "Angiogenic tip",
         X_escape = "X-escape (positive control)", Stress = "Dissociation stress",
         peri_cap_ratio = "Pericytes per capillary EC", cap_share = "Capillary share of ECs", EC_share = "EC share of cells")

# ---- Figure 2 ----
m1 <- A1$main |> filter(program %in% c(PRIMARY, "X_escape")) |>
  mutate(label = factor(lab[program], levels = rev(lab[c(PRIMARY, "X_escape")])))
d1 <- A1$per_ds |> filter(program %in% c(PRIMARY, "X_escape")) |> left_join(dsn, by = "dataset_id") |>
  mutate(label = factor(lab[program], levels = levels(m1$label)))
f2 <- ggplot() +
  annotate("rect", xmin = -MARGIN, xmax = MARGIN, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(data = d1, aes(est, label, size = n_w + n_m), shape = 1, colour = "grey55",
             position = position_jitter(width = 0, height = 0.15, seed = 1)) +
  geom_errorbar(data = m1, aes(xmin = lo, xmax = hi, y = label), width = 0, orientation = "y", linewidth = 0.6) +
  geom_errorbar(data = m1, aes(xmin = lo90, xmax = hi90, y = label), width = 0.25, orientation = "y", linewidth = 0.3) +
  geom_point(data = m1, aes(est, label), shape = 18, size = 3.2) +
  scale_size_area(max_size = 3, name = "Donors") +
  labs(x = "Women minus men (donor-level SD units)", y = NULL,
       caption = "Diamonds: pooled estimate with 95% CI (thick) and 90% CI (thin). Grey band: equivalence bounds (±0.8 SD). Circles: single datasets.")
save_fig(f2, "Figure2_capillary_programs", 95)

# ---- Figure 3 ----
f3a <- ggplot(A1$perm |> filter(program %in% PRIMARY) |> mutate(label = lab[program]), aes(p)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey70", colour = "white") +
  facet_wrap(~ label, nrow = 1) + labs(x = "Meta-analysis p value, shuffled sex labels", y = "Permutations")
g <- A1$gmeta; if (!"chr" %in% names(g)) g$chr <- NA
g <- g |> mutate(cls = case_when(chr == "X" ~ "X", chr == "Y" ~ "Y", TRUE ~ "Autosome"))
f3b <- ggplot(g, aes(logFC, -log10(p), colour = cls)) + geom_point(size = 0.5, alpha = 0.6) +
  geom_text(data = g |> filter(FDR < 0.05) |> arrange(p) |> head(15), aes(label = gene), size = 2.2,
            vjust = -0.6, show.legend = FALSE) +
  scale_colour_manual(values = c(Autosome = "grey60", X = W, Y = M), name = NULL) +
  labs(x = "Women minus men (log2, pooled)", y = expression(-log[10]~p))
save_fig(patchwork::wrap_plots(f3a, f3b, ncol = 1, heights = c(1, 1.6)) + patchwork::plot_annotation(tag_levels = "A"),
         "Figure3_calibration_genomewide", 150)

# ---- Figure 4 ----
m2 <- bind_rows(A2$comp |> mutate(fam = "Composition"),
                A2$lr |> mutate(fam = "Capillary-pericyte signaling")) |>
  mutate(label = ifelse(program %in% names(lab), lab[program], program))
f4 <- ggplot(m2, aes(est, reorder(label, est))) +
  annotate("rect", xmin = -MARGIN, xmax = MARGIN, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.6) +
  geom_point(shape = 18, size = 3) + facet_grid(fam ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Women minus men (donor-level SD units)", y = NULL)
save_fig(f4, "Figure4_structure_signaling", 90)

# ---- Figure 5 ----
h <- A4$by_group |> filter(measure %in% c("X_only", "X_plus_Y")) |>
  mutate(measure = recode(measure, X_only = "X copy alone", X_plus_Y = "X + Y combined"),
         star = ifelse(FDR < 0.05, "*", ""))
f5 <- ggplot(h, aes(group, pair, fill = est)) + geom_tile(colour = "white") +
  geom_text(aes(label = paste0(sprintf("%+.2f", est), star)), size = 2.3) +
  facet_wrap(~ measure) +
  scale_fill_gradient2(low = M, mid = "white", high = W, midpoint = 0, name = "Women - men\n(log2)") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_fig(f5, "Figure5_XY_paralogs", 90)

# ---- Tables (CSV) ----
dir.create("tables", showWarnings = FALSE)
write.csv(A1$main, "tables/Table2_aim1_programs.csv", row.names = FALSE)
write.csv(A1$per_ds |> left_join(dsn, by = "dataset_id"), "tables/TableS_aim1_per_dataset.csv", row.names = FALSE)
write.csv(A1$gmeta |> arrange(p), "tables/TableS_genomewide_meta.csv", row.names = FALSE)
write.csv(bind_rows(A2$comp, A2$lr), "tables/Table3_aim2.csv", row.names = FALSE)
write.csv(bind_rows(A3$main, A3$sens), "tables/TableS_aim3_age.csv", row.names = FALSE)
write.csv(A4$by_group, "tables/Table4_aim4_by_group.csv", row.names = FALSE)
write.csv(A4$spec, "tables/TableS_aim4_specificity.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
message("Figures in ", file.path(getwd(), "figures"), "; tables in ", file.path(getwd(), "tables"))
