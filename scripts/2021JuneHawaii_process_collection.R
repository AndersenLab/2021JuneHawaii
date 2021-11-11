library(tidyverse)
library(easyfulcrum)
#devtools::install_github("AndersenLab/easyfulcrum")

# define project directory
dir = "~/repos/2021JuneHawaii"

# make dir structure for processing
makeDirStructure(startdir = dir)

# read fulcrum data
raw_fulc <- readFulcrum(dir = dir)

# process the data
proc_fulc <- procFulcrum(data = raw_fulc)

# check temps, only four measures with same ambient humidity reading, they should not be a problem.
# no need to run fixTemperatures()
flag_temp <- checkTemperatures(data = proc_fulc, return_flags = TRUE)

# join the fulcrum data
join_fulc <- joinFulcrum(data = proc_fulc, select_vars = TRUE) %>%
  dplyr::filter(!(s_label %in% c("S-14574", "S-14444"))) %>% # remove records with duplicated S-label S-14574 and mis genotyped S-label
  dplyr::mutate(substrate_other = as.character(substrate_other),  # fix data classes
                gps_speed = as.numeric(gps_speed),  # fix data classes
                gps_vertical_accuracy = as.numeric(gps_vertical_accuracy)) # fix data classes

# checking the join
flag_join <- checkJoin(data = join_fulc, return_flags = TRUE) 
# S-14574 is duplicated in Fulcrum, this s-label was added by Erik twice.
#   isolation: e96be254-04b0-47b7-8381-3ecbcf854064, c_label: 9e05fa86-74b8-4795-acf1-332d7e844460 C-6535
#   isolation: c6e5f396-4797-4f3b-93ce-12aac0fd9df5, c_label: 041c1215-4877-476b-8d09-8f33c794427b C-6536
# 2021-11-10: My guess is Erik may have accidentally scanned S-14574 instead of S-14573 for C-6536,
#   b/c S-14573 is missing and numeric sequences for other s_plates and C-plates match that expectation.
#   Safest solution is to remove both instances of S-14574. My theory on S-14573 seems to be validated b/c
#   the checkGenotypes function found S-14573 in the genotyping data.

# annotate the collections
anno_fulc <- annotateFulcrum(data = join_fulc, dir = NULL, select_vars = TRUE)

# Read genotyping sheet
raw_geno_nema <- readGenotypes(gsKey = c("1DpsOa6C6zMx9qxq8BthusiiVK8Lw-7S5bKTdcMoaN7o"),
                               col_types = "cDDdcdcddddddDcDDdcdcdddddddcdcccddccc") %>%
  dplyr::filter(!(s_label %in% c("S-14573", "S-14574", "S-14444"))) # fix S-laebls with errors

# process the genotyping sheet
proc_geno_nema <- checkGenotypes(geno_data = raw_geno_nema, fulc_data = anno_fulc, 
                                 return_geno = TRUE, return_flags = FALSE, profile = "nematode") 
# S-14573 is in genotyping sheet but not fulcrum. This is prbably b/c Erik accidentally scanned S-14574 in 
#   its place. The isolate has an SSU band but not an ITS2 band I will drop S-14573 and S-14574.
# S-14444 is lost, Robyn said it failed to sequence but could not be re-lysed. We are dropping from data.

# join genotype data with Fulcrum data
join_genofulc_nema <- joinGenoFulc(geno = proc_geno_nema, fulc = anno_fulc, dir = NULL, select_vars = TRUE)

# Process photos
final_data_nema <- procPhotos2(dir = dir, data = join_genofulc_nema,
                              max_dim = 500, overwrite = TRUE,
                              pub_url = "https://storage.googleapis.com/elegansvariation.org/photos/isolation/fulcrum/",
                              CeaNDR = TRUE)

# make the species sheet for CeNDR
sp_sheet <- makeSpSheet(data = final_data_nema, dir = dir)
# lots of issues to fix here, what's up with landscapes?

# Make final report
generateReport(data = final_data_nema, dir = dir, profile = "nematode")
