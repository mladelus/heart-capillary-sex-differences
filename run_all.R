# run_all.R - run every step in order. From the repository folder:  source("run_all.R")
# Steps that download data skip donors already saved, so the script can be re-run after an interruption.
# Set the data folder first if you do not want the default (~/Documents/heart_capillary_project):
#   Sys.setenv(HCAP_DATA = "/path/to/folder")
repo <- normalizePath(getwd())
run  <- function(f) { message("\n=========== ", f, " ==========="); source(file.path(repo, "R", f), echo = FALSE) }
for (f in c("00_setup.R", "01_inventory.R", "02_download.R", "03_ec_subtypes_pseudobulk.R",
            "04_aim1_capillary_programs.R", "05_aim2_structure.R", "06_aim3_age.R",
            "07_aim4_sex_chromosomes.R", "08_figures.R")) run(f)
