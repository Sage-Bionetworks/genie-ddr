library(purrr)
library(here)
library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

samp_aug_fs <- readr::read_rds(
  here('data', 'combined', 'samp_aug_first_sample.rds')
)


# We'll calculate this as a vector then put it back in.
# rowSums() is just a huge pain with pipes, but it's fast.
any_ddr <- samp_aug_fs %>%
  select(any_of(custom_ddr_list())) %>%
  rowSums(., na.rm = T) %>%
  magrittr::is_greater_than(., 0)
samp_aug_fs %<>% mutate(any_ddr = any_ddr)

tumor_types <- mskcc.oncotree::get_tumor_types()
samp_aug_fs %<>%
  left_join(
    .,
    select(tumor_types, oncotree_code, tissue),
    by = "oncotree_code"
  )

samp_aug_fs %<>%
  mutate(
    custom_cancer_type = case_when(
      oncotree_code %in% c("UNKNOWN", "CUP", "CUPNOS") ~ "Unknown",
      # Three types for bladder:
      oncotree_code %in% "UTUC" ~ "Bladder (UTUC)",
      oncotree_code %in% "BLCA" ~ "Bladder (BLCA)",
      tissue %in% "Bladder/Urinary Tract" ~ "Bladder (other)",
      # For some reason glioma's dont have a tissue type...
      # My fix is to use the cancer type for those, as much as I hate that.
      is.na(tissue) ~ cancer_type,
      TRUE ~ tissue
    )
  )

# We're losing some important groups doing this, but I want to see the plot:
samp_aug_fs %<>%
  group_by(custom_cancer_type) %>%
  mutate(n_custom_cancer_type = n()) %>%
  ungroup(.) %>%
  filter(n_custom_cancer_type >= 10)

# Make a feature for cancer type with an n:
samp_aug_fs %<>%
  group_by(custom_cancer_type) %>%
  mutate(n_ct = n()) %>%
  ungroup(.) %>%
  # putting these in order with weird sorting so I don't have to mess with the n's later on:
  arrange(
    desc(str_detect(custom_cancer_type, "BLCA")),
    desc(str_detect(custom_cancer_type, "UTUC")),
    desc(str_detect(custom_cancer_type, "Bladder \\(other\\)")),
    custom_cancer_type
  )

samp_aug_fs %<>%
  mutate(
    # leaving the paren open here on purpose:
    custom_cancer_type = glue('{custom_cancer_type} ({n_ct}'),
    custom_cancer_type = fct_inorder(custom_cancer_type),
    custom_cancer_type = fct_rev(custom_cancer_type)
  )

cts <- samp_aug_fs %>%
  group_by(custom_cancer_type) %>%
  rename(ddr = any_ddr) %>%
  summarize(
    across(
      .cols = c(any_of(custom_ddr_list()), ddr),
      .fns = list(
        # proportion of samples positive out of those tested:
        pos = \(x) {
          sum(x, na.rm = T) / sum(!is.na(x))
        },
        # proportion of samples tested:
        test = \(x) {
          sum(!is.na(x)) / n()
        }
      )
    )
  )

# This block turns the ddr features into text and then gets rid of them:
cts %<>%
  mutate(
    custom_cancer_type_pos_lab = glue(
      "{custom_cancer_type}, {formatC(ddr_pos*100, digits = 0, format = 'f')}%)"
    ),
    custom_cancer_type_test_lab = glue(
      "{custom_cancer_type}, {formatC(ddr_test*100, digits = 0, format = 'f')}%)"
    ),
    custom_cancer_type_pos_lab = fct_inorder(custom_cancer_type_pos_lab),
    custom_cancer_type_test_lab = fct_inorder(custom_cancer_type_test_lab)
  )


cts_long <- cts %>%
  pivot_longer(
    cols = -contains("custom_cancer_type"),
    names_to = "feature",
  )

cts_long %<>%
  separate_wider_delim(
    feature,
    delim = "_",
    names = c('feature', 'metric'),
    too_many = "error",
    too_few = 'error'
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
    desc(feature %in% 'any_DDR'),
    feature
  ) %>%
  mutate(feature = fct_inorder(feature))

gg_pos <- plot_ddr_heatmap(
  filter(cts_long, feature != "any_DDR"),
  plot_title = "Proportion of first samples altered by cancer type",
  plot_subtitle = "Denominator is the number of samples TESTED for individual genes"
)
ggsave(
  gg_pos,
  height = 15,
  width = 10,
  filename = here('output', 'fig', '1A_fig_bladder_manu_main_genie_ddr_pos.pdf')
)


gg_test <- ggplot(
  (cts_long %>% filter(!(feature %in% "any_DDR"))),
  aes(x = feature, y = custom_cancer_type_test_lab, fill = test)
) +
  geom_tile(color = 'gray20') +
  scale_fill_viridis_c(
    option = "viridis",
    begin = 0.2,
    end = 0.8,
    name = "Proportion tested"
  ) +
  scale_y_discrete(position = 'right', expand = c(0, 0)) +
  scale_x_discrete(position = 'top', expand = c(0, 0)) +
  theme_bw() +
  labs(
    title = "Proportion of first samples tested for each cancer type",
    subtitle = "Denominator is the number of samples (first sample for each person)"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = 'top'
  )

ggsave(
  gg_test,
  height = 15,
  width = 10,
  filename = here(
    'output',
    'fig',
    'S1_fig_bladder_manu_main_genie_ddr_test.pdf'
  )
)
