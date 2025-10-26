# Description:  Assess the impact of oncoKB filtering, save files with
#  only the variants which pass an oncoKB pass.

library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

# here onco just refers to it having the oncogenic tag, not being limited to those alterations.
mut_onco <- readr::read_tsv(
  here('data', 'genomic', 'mut_onco.txt'),
  show_col_types = F
)

samp_aug_proto <- readr::read_rds(
  here('data', 'combined', 'samp_aug_proto.rds')
)


# Just some print outs for the analyst - hoping for near 100% in all:
genomic_anno_msg_help(mut_onco)

lev_onco <- c(
  "Oncogenic",
  "Likely Oncogenic",
  "Resistance",
  "Likely Neutral",
  "Inconclusive",
  "Unknown"
)

mut_onco %<>%
  rename_all(tolower) %>%
  mutate(
    mut_vaf = t_alt_count / (t_alt_count + t_ref_count),
    oncogenic = factor(oncogenic, levels = lev_onco)
  )

# separate script or something probably best:
#
# readr::write_tsv(
#   mut_onco,
#   here('data', 'genomic', 'mut_onco_processed.rds')
# )

onco_gene_feat <- mut_onco %>%
  filter(oncogenic %in% c('Likely Oncogenic', 'Oncogenic')) %>%
  select(
    sample_id = tumor_sample_barcode,
    hugo_symbol
  ) %>%
  group_by(sample_id, hugo_symbol) %>%
  summarize(altered = T, .groups = 'drop')

# This object is larger, but we'll be able to save memory later on:
onco_gene_feat_mat <- long_dat_to_mat_helper(onco_gene_feat)
onco_gene_feat_mat[is.na(onco_gene_feat_mat)] <- FALSE


testing_mat <- readr::read_rds(
  here('data', 'genomic', 'matrix_testing_by_hugo.rds')
)

onco_comb_mat <- combine_testing_and_alteration_matrices(
  testing_mat,
  onco_gene_feat_mat
)

rm(testing_mat)
rm(onco_gene_feat_mat)
gc()

onco_gene_feat <- as_tibble(onco_comb_mat, rownames = "sample_id")


cli_abort("Works up to here with no memory issues.")

# here's a thought:  the drugs are listed under level of evidence.  could scan for cisplatin there.
any_gene_feat <- mut_onco %>%
  # No oncogenic filter at all here
  select(
    sample_id = tumor_sample_barcode,
    hugo_symbol
  ) %>%
  group_by(sample_id, hugo_symbol) %>%
  summarize(altered = T, .groups = 'drop')

rm(mut_onco) # large memory use - don't need it now.

any_gene_feat_mat <- long_dat_to_mat_helper(any_gene_feat)
any_gene_feat_mat[is.na(any_gene_feat_mat)] <- FALSE


onco_gene_feat <- onco_gene_feat %>%
  select(sort(tidyselect::peek_vars())) %>%
  select(sample_id, everything())

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


cli_abort("stop")
# Repeat those steps for the any_gene_feat version:

any_gene_feat %<>%
  full_join(
    samp_skel,
    .,
    by = c('sample_id', 'hugo_symbol')
  )

any_gene_feat %<>%
  replace_na(list(tested = F, altered = F))
any_gene_feat %<>%
  mutate(
    feature = case_when(
      altered ~ T, # including cases when untested, at least for now.
      tested ~ F,
      T ~ NA
    )
  )
any_gene_feat <- any_gene_feat %>%
  select(sample_id, hugo_symbol, feature) %>%
  pivot_wider(
    names_from = 'hugo_symbol',
    values_from = feature
  )
any_gene_feat <- any_gene_feat %>%
  select(sort(tidyselect::peek_vars())) %>%
  select(sample_id, everything())


samp_aug_onco <- left_join(
  samp_aug_proto,
  onco_gene_feat,
  by = 'sample_id'
)

samp_aug_any <- left_join(
  samp_aug_proto,
  any_gene_feat,
  by = 'sample_id'
)


readr::write_rds(
  samp_aug_onco,
  here('data', 'combined', 'samp_aug_onco.rds')
)

readr::write_rds(
  samp_aug_onco,
  here('data', 'combined', 'samp_aug_any.rds')
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
