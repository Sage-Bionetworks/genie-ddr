samp_aug_fs <- readr::read_rds(
  here('data', 'combined', 'samp_aug_first_sample.rds')
)

custom_ddr_list <- c(
  'ERCC2', 'ERCC5', 
  'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
  'ATR', 'FANCC'
) %>% sort

samp_aug_fs %>%
#  select(any_of(custom_ddr_list)) %>%
  purrr::map(
    .x = .,
    .f = ~sum(is.na(.x))
  )

samp_aug_fs %<>% 
  mutate(
    any_ddr = rowSums(
      mutate(
        select(., any_of(custom_ddr_list),
               across(.cols = everything(), .fns = ~replace_na(.x, F)))
      )
    )
  )

samp_aug_fs %>% count(cancer_type)

# We're losing some important groups doing this, but I want to see the plot:

samp_aug_fs %<>%
  group_by(cancer_type) %>%
  mutate(n_cancer_type = n()) %>%
  ungroup(.) %>%
  filter(n_cancer_type > 99)

cts <- samp_aug_fs %>%
  group_by(cancer_type) %>%
  summarize(
    across(
      # .cols = c(any_of(custom_ddr_list), any_ddr),
      .cols = c(any_of(custom_ddr_list)),
      .fns = \(x) sum(x) / n()
    )
  )

cts_long <- cts %>%
  pivot_longer(
    cols = -cancer_type,
    names_to = "feature",
    values_to = "prop_ddr"
  ) %>%
  mutate(cancer_type = factor(cancer_type),
         cancer_type = fct_rev(cancer_type)
  )

gg <- ggplot(
  cts_long,
  aes(x = feature, y = cancer_type, fill = prop_ddr)
) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  scale_y_discrete(position = 'right') +
  scale_x_discrete(position = 'top') + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0)
  )

ggsave(
  gg, height = 10, width = 5,
  filename = here('output', 'fig', 'ddr_alt_cancer_type.pdf')
)

