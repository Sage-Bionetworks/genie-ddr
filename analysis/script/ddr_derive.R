
library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

samp_aug_fs <- readr::read_rds(
  here('data', 'combined', 'samp_aug_first_sample.rds')
)

custom_ddr_list <- c(
  'ERCC2', 'ERCC5', 
  'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
  'ATR', 'FANCC'
) %>% sort

# Can check how often these were tested easily now:
# samp_aug_fs %>%
#   select(any_of(custom_ddr_list)) %>%
#   purrr::map(
#     .x = .,
#     .f = ~sum(is.na(.x))
#   )

# We'll calculate this as a vector then put it back in.
# rowSums() is just a huge pain with pipes, but it's fast.
any_ddr <- samp_aug_fs %>%
  select(any_of(custom_ddr_list)) %>%
  rowSums(., na.rm = T) %>%
  magrittr::is_greater_than(.,0)
samp_aug_fs %<>% mutate(any_ddr = any_ddr)

samp_aug_fs %>% count(cancer_type)
samp_aug_fs %>% filter(cancer_type %in% "Bladder Cancer") %>% count(oncotree_code)

# We're losing some important groups doing this, but I want to see the plot:
samp_aug_fs %<>%
  group_by(cancer_type) %>%
  mutate(n_cancer_type = n()) %>%
  ungroup(.) %>%
  filter(n_cancer_type >= 10)

# Make a feature for cancer type with an n:
samp_aug_fs %<>%
  group_by(cancer_type) %>%
  mutate(n_ct = n()) %>%
  ungroup(.) %>%
  mutate(
    cancer_type = glue('{cancer_type} ({n_ct})'),
    cancer_type = factor(cancer_type),
    cancer_type = fct_rev(cancer_type)
  ) 

cts <- samp_aug_fs %>%
  group_by(cancer_type) %>%
  rename(ddr = any_ddr) %>%
  summarize(
    across(
      .cols = c(any_of(custom_ddr_list), ddr),
      # .cols = c(any_of(custom_ddr_list)),
      .fns = list(
        # proportion of samples positive out of those tested:
        pos = \(x) {sum(x, na.rm = T) / sum(!is.na(x))},
        # proportion of samples tested:
        test = \(x) {sum(!is.na(x)) / n()}
      )
    )
  )

cts_long <- cts %>%
  pivot_longer(
    cols = -cancer_type,
    names_to = "feature",
  ) 

cts_long %<>%
  separate_wider_delim(
    feature,
    delim = "_",
    names = c('feature', 'metric'),
    too_many = "error", too_few = 'error'
  ) %>%
  pivot_wider(
    names_from = 'metric',
    values_from = 'value'
  ) 

cts_long %<>%
  mutate(
    feature = if_else(feature %in% 'ddr', 'any_DDR', feature)
  ) %>%
  arrange(
    desc(feature %in% 'any_DDR'), feature
  ) %>%
  mutate(feature = fct_inorder(feature))
    

gg_pos <- ggplot(
  cts_long,
  aes(x = feature, y = cancer_type, fill = pos)
) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma", name = "Proportion altered") + 
  scale_y_discrete(position = 'right') +
  scale_x_discrete(position = 'top') + 
  theme_bw() + 
  labs(
    title = "Proportion of first samples altered by cancer type",
    subtitle = "Denominator is the number of samples TESTED for individual genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0),
    legend.position = 'top'
  )

ggsave(
  gg_pos, height = 15, width = 10,
  filename = here('output', 'fig', 'ddr_alt_cancer_type.pdf')
)





gg_test <- ggplot(
  (cts_long %>% filter(!(feature %in% "any_DDR"))),
  aes(x = feature, y = cancer_type, fill = test)
) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "viridis", begin = 0.2, end = 0.8, name = "Proportion tested") + 
  scale_y_discrete(position = 'right') +
  scale_x_discrete(position = 'top') + 
  theme_bw() + 
  labs(
    title = "Proportion of first samples tested for each cancer type",
    subtitle = "Denominator is the number of samples (first sample for each person)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0),
    legend.position = 'top'
  )

ggsave(
  gg_test, height = 15, width = 10,
  filename = here('output', 'fig', 'ddr_test_cancer_type.pdf')
)
