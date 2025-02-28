# Description:  Assess the impact of oncoKB filtering, save files with
#  only the variants which pass an oncoKB pass.

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

mut_onco <- readr::read_tsv(
  here('data', 'genomic', 'mut_onco.txt'),
  show_col_types = F
)
cna_onco <- readr::read_tsv(
  here('data', 'genomic', 'cna_onco.txt'),
  show_col_types = F
)
sv_onco <- readr::read_tsv(
  here('data', 'genomic', 'fus_onco.txt'),
  show_col_types = F
)

# For CNAs we focus only on the highly amplified cases, taking the decisions from the breast manuscript as a good starting point if nothing else.
dft_cna_raw <- readr::read_tsv(
  here('data', 'genomic', 'cna_ready_to_annotate.txt'),
  show_col_types = F
) 
cna_onco <- bind_cols(
  # The sample ID that came back is garbage.
  select(cna_onco, -SAMPLE_ID),
  select(dft_cna_raw,  SAMPLE_ID = Tumor_Sample_Barcode)
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
anno_msg_help(cna_onco)
anno_msg_help(sv_onco)
