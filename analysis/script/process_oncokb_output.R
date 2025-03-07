# Description:  Assess the impact of oncoKB filtering, save files with
#  only the variants which pass an oncoKB pass.

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

mut_onco <- readr::read_tsv(
  here('data', 'genomic', 'mut_onco.txt'),
  show_col_types = F
)

samp <- readr::read_tsv(
  file = here('data-raw', 'main_genie', 'data_clinical_sample.txt'),
  comment = "#"
)

clin <- readr::read_tsv(
    file = here('data-raw', 'main_genie', 'data_clinical_patient.txt'),
    comment = "#"
  )
  
  


anno_msg_help <- function(dat) {
  dat_name <- deparse(substitute(dat))
  
  tab <- tabyl(dat, ANNOTATED) %>%
    filter(ANNOTATED) # T/F so this grabs the annotated pct line.
  
  n_levels_oncogenic <- count(dat, ONCOGENIC) %>%
    nrow(.)
  
  cli::cli_alert_info(
    text = glue("{dat_name}: Of the {tab$n} rows, {round(tab$percent*100,1)}% were annotated.")
  )
  
  # just a check that everything wasn't tagged the same way,
  #   or something equally nonsensical.
  cli::cli_alert_info(
    text = glue("{dat_name}: {n_levels_oncogenic} levels for ONCOGENIC (just hoping >= 2).")
  )
}

# Just some print outs for the analyst - hoping for near 100% in all:
anno_msg_help(mut_onco)


lev_onco <- c("Oncogenic", "Likely Oncogenic",
              "Resistance",
              "Likely Neutral",
              "Inconclusive", "Unknown")

mut_onco %<>%
  rename_all(tolower) %>%
  mutate(
    mut_vaf = t_alt_count / (t_alt_count + t_ref_count),
    oncogenic = factor(oncogenic, levels = lev_onco)
  )

readr::write_tsv(
  mut_onco,
  here('data', 'genomic', 'mut_onco_processed.rds')
)


# Because I have yet to see resistance:
mut_onco %>% count(oncogenic)
mut_onco %>% 
  filter(oncogenic %in% "Resistance") %>%
  count(hugo_symbol) %>%
  print(n = 500)
# not that many genes, none on the core DDR list, so moving on.

# here's a thought:  the drugs are listed under level of evidence.  could scan for cisplatin there.
onco_gene_feat <- mut_onco %>%
  filter(oncogenic %in% c('Likely Oncogenic', 'Oncogenic')) %>%
  select(
    sample_id = tumor_sample_barcode,
    hugo_symbol
  ) %>%
  group_by(sample_id, hugo_symbol) %>%
  summarize(altered = T, .groups = 'drop')

readr::write_rds(
  onco_gene_feat,
  here('data', 'genomic', 'onco_mut_long.rds')
)


gi <- read_tsv(
  here('data-raw', 'main_genie', 'genomic_information.txt')
)
# Based on the best information we have now (inadequate) we believe includeInPanel
#.  is true when the region would have been reported if positive for mutations.
# includeInPanel = F is for regions which were in the vendor's file, or part of 
#   a copy number probe, or something like that.
# Also note that we have lots of introns and intergenic regions in the F parts:
# gi %>% janitor::tabyl(includeInPanel, Feature_Type)

gi_simp <- gi %>%
  filter(includeInPanel) %>%
  group_by(SEQ_ASSAY_ID, Hugo_Symbol) %>%
  summarize(tested = T, .groups = 'drop') %>%
  rename_all(tolower)
  
samp_skel <- samp %>%
  rename_all(tolower) %>%
  select(sample_id, seq_assay_id) 

if (length(setdiff(unique(samp_skel$seq_assay_id), gi_simp$seq_assay_id)) > 0) {
  cli::cli_alert_warning("Some assays. in the sample file were not found in the genomic infomation (BAM) file.  This will result in these counting as 'not tested' for a gene.  Alert data team if you see this.")
}


# Joins are probably a memory-inefficient way to do this, but I don't care, it works for now.
samp_skel <- left_join(samp_skel, gi_simp, by = 'seq_assay_id',
                       relationship = 'many-to-many')

onco_gene_feat %<>%
  full_join(
    samp_skel,
    .,
    by = c("sample_id", 'hugo_symbol')
  ) 

rm(samp_skel) # big file.

onco_gene_feat %<>%
  replace_na(list(tested = F, altered = F))

untested_but_altered <- onco_gene_feat %>%
  filter(!tested & altered)

if (nrow(untested_but_altered) > 0) {
  n_untest_alt <- nrow(untested_but_altered)
  prop_untest_alt <- n_untest_alt / nrow(onco_gene_feat)
  cli::cli_alert_warning("{n_untest_alt} ({round(prop_untest_alt)}% of rows) alterations detected which were not tested according to genomic information file.")
  readr::write_csv(
    # seq assay id is all na because we joined the assay to the testing data (not the altered data).
    select(untested_but_altered, -seq_assay_id),
    here('analysis', 'explore', 'untested_altered_variants.csv')
  )
}

onco_gene_feat %<>%
  mutate(
    feature = case_when(
      altered ~ T, # including cases when untested, at least for now.
      tested ~ F,
      T ~ NA
    )
  )

onco_gene_feat %<>%
  select(sample_id, hugo_symbol, feature) %>%
  pivot_wider(
    names_from = 'hugo_symbol',
    values_from = feature
  )

onco_gene_feat <- onco_gene_feat %>%
  select(sort(tidyselect::peek_vars())) %>%
  select(sample_id, everything())





# While we're here we'll create a version with all the clinical and genomic features:
samp %<>% rename_all(tolower)
clin %<>% rename_all(tolower)

samp %<>%
  select(-age_at_seq_report_days) %>%
  select(
    sample_id, patient_id, 
    age_at_seq_report:sample_class
  ) 

clin %<>% 
  # secondary/tertiary race are useless,so:
  select(patient_id, sex, primary_race, ethnicity:year_death)

# will need to fix up some of the age variables, but it's not important now so we'll wait.

samp_aug <- left_join(
  samp, clin, by = 'patient_id'
)










samp_aug <- left_join(
  samp_aug, onco_gene_feat, by = 'sample_id'
)

samp_aug %<>%
  mutate(
    age_at_seq_report_mod = case_when(
      age_at_seq_report %in% "<18" ~ 17,
      # so annoying, lots of these people are probably 89 not 90.
      age_at_seq_report %in% ">89" ~ 89,
      age_at_seq_report %in% "Unknown" ~ NA_real_,
      T ~ as.numeric(age_at_seq_report)
    )
  ) %>%
  relocate(
    age_at_seq_report_mod, .after = age_at_seq_report
  )

# This genomic dataset is positive/else only - no difference between NA and False for now
samp_aug %<>%
  mutate(
    across(
      .cols = -c(sample_id:year_death),
      .fns = ~replace_na(.x, replace = F)
    )
  ) 

readr::write_rds(
  samp_aug,
  here('data', 'combined', 'samp_aug.rds')
)


samp_aug_first_sample <- samp_aug %>%
  group_by(patient_id) %>%
  # if there are ties for both year and age we'll just 
  #   pick arbitrarily (sample_id).
  arrange(seq_year, age_at_seq_report_mod, sample_id) %>%
  slice(1) %>%
  ungroup(.)

readr::write_rds(
  samp_aug_first_sample,
  here('data', 'combined', 'samp_aug_first_sample.rds')
)
  




