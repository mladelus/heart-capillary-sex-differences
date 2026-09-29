install.packages(c("metafor", "patchwork"))
BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)

setwd("~/Desktop/heart-capillary-sex-differences")
source("~/Desktop/heart-capillary-sex-differences/R/00_setup.R")
source("~/Desktop/heart-capillary-sex-differences/R/01_inventory.R")
