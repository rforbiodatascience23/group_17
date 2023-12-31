

oad data:
  
  {r}

protein_data <- read.table("../data/count_matrix_1_77.tsv", header = T) |>
  pivot_longer(cols = -id, names_to = 'variable', values_to = 'value') |> 
  pivot_wider(names_from = id, values_from = value) |>
  rename(id = variable)

clinical_data <- read.table("../data/clinical_data.tsv", sep = "\t", header = T)

Wrangling

{r}

full_data <- clinical_data |>
  mutate(
    DGEA_HER2 = case_when(
      pam50 == "HER2-enriched" ~ 1, 
      TRUE ~ 0
    ), 
    DGEA_Luminal_A = case_when(
      pam50 == "Luminal A" ~ 1, 
      TRUE ~ 0
    ), 
    DGEA_Luminal_B = case_when(
      pam50 == "Luminal B" ~ 1, 
      TRUE ~ 0
    ),
    DGEA_Basal_like = case_when(
      pam50 == "Basal-like" ~ 1, 
      TRUE ~ 0)) |>
  left_join(protein_data, by = "id")


#HER2 positive
HER2 <- full_data |>
  select(c(DGEA_HER2, starts_with("NP"))) |>
  pivot_longer(cols = starts_with("NP"),
               names_to = "protein",
               values_to = "value") |>
  group_by(protein) |>
  nest() |> 
  mutate(LM = map(.x = data,
                  .f = ~lm(formula = DGEA_HER2 ~ value,
                           data = .x))) |>
  mutate(LM_tidy = map(.x = LM,
                       .f = ~broom::tidy(x = .x,
                                         conf.int = TRUE,
                                         conf.level = 0.95))) |>
  unnest(LM_tidy) |>
  filter(term == "value") |>
  select(-c(data, LM, term)) |>
  mutate(adjusted_BH = p.adjust(p.value, n = 12554, method = "BH"), 
         adjusted_fdr = p.adjust(p.value, n = 12554, method = "fdr"), 
         log2 = -log2(p.value)) |>
  arrange(p.value)


#Luminal A
Luminal_A <- full_data |>
  select(c(DGEA_Luminal_A, starts_with("NP"))) |>
  pivot_longer(cols = starts_with("NP"),
               names_to = "protein",
               values_to = "value") |>
  group_by(protein) |>
  nest() |> 
  mutate(LM = map(.x = data,
                  .f = ~lm(formula = DGEA_Luminal_A ~ value,
                           data = .x))) |>
  mutate(LM_tidy = map(.x = LM,
                       .f = ~broom::tidy(x = .x,
                                         conf.int = TRUE,
                                         conf.level = 0.95))) |>
  unnest(LM_tidy) |>
  filter(term == "value") |>
  select(-c(data, LM, term)) |>
  mutate(adjusted_BH = p.adjust(p.value, n = 12554, method = "BH"), 
         adjusted_fdr = p.adjust(p.value, n = 12554, method = "fdr")) |>
  arrange(p.value)

#Luminal B
Luminal_B <- full_data |>
  select(c(DGEA_Luminal_B, starts_with("NP"))) |>
  pivot_longer(cols = starts_with("NP"),
               names_to = "protein",
               values_to = "value") |>
  group_by(protein) |>
  nest() |> 
  mutate(LM = map(.x = data,
                  .f = ~lm(formula = DGEA_Luminal_B ~ value,
                           data = .x))) |>
  mutate(LM_tidy = map(.x = LM,
                       .f = ~broom::tidy(x = .x,
                                         conf.int = TRUE,
                                         conf.level = 0.95))) |>
  unnest(LM_tidy) |>
  filter(term == "value") |>
  select(-c(data, LM, term)) |>
  mutate(adjusted_BH = p.adjust(p.value, n = 12554, method = "BH"), 
         adjusted_fdr = p.adjust(p.value, n = 12554, method = "fdr")) |>
  arrange(p.value)


#Basal like
Basal_like <- full_data |>
  select(c(DGEA_Basal_like, starts_with("NP"))) |>
  pivot_longer(cols = starts_with("NP"),
               names_to = "protein",
               values_to = "value") |>
  group_by(protein) |>
  nest() |> 
  mutate(LM = map(.x = data,
                  .f = ~lm(formula = DGEA_Basal_like ~ value,
                           data = .x))) |>
  mutate(LM_tidy = map(.x = LM,
                       .f = ~broom::tidy(x = .x,
                                         conf.int = TRUE,
                                         conf.level = 0.95))) |>
  unnest(LM_tidy) |>
  filter(term == "value") |>
  select(-c(data, LM, term)) |>
  mutate(adjusted_BH = p.adjust(p.value, n = 12554, method = "BH"), 
         adjusted_fdr = p.adjust(p.value, n = 12554, method = "fdr")) |>
  arrange(p.value)

Save:
  
  {r}
write_tsv(HER2, "../data/diff_exp_HER2.tsv")
write_tsv(Luminal_A, "../data/diff_exp_Luminal_A.tsv")
write_tsv(Luminal_B, "../data/diff_exp_Luminal_B.tsv")
write_tsv(Basal_like, "../data/diff_exp_Basal_like.tsv")

For Carina, a template for volcano plots:
  
  Tip: Start from the top. The first lines are the most important ones. Afterwards it is just krymmel.

{r}

HER2 <- read.table()

HER2 <- HER2 |>
  mutate(color_col = case_when(
    adjusted_fdr < 0.05 & estimate > 0 ~ "Significant positive influence", 
    adjusted_fdr < 0.05 & estimate < 0~ "Significant negative influence",
    TRUE ~"Not significant"
  ))

library(ggrepel)

plot <- ggplot(HER2, aes(x = estimate, y = log2)) +
  geom_point(aes(color = color_col)) +
  labs(title = "HER2", 
       x = "Estimate", 
       y = "-Log2(FDR-adjusted P-value)") +
  geom_hline(yintercept = 18, linetype = "dashed", color = "black" ) +
  geom_hline(yintercept = (2^-0.8), linetype = "dashed", color = "black" ) +
  scale_color_manual(values = c("Significant positive influence" = "orange", 
                                "Significant negative influence" = "#69b3a2",
                                "No significant association" = "grey")) +
  geom_text_repel(data = HER2 |>filter(adjusted_fdr < 0.05),
                  aes(label = protein), hjust = 0, vjust = 0, size = 3, 
                  color = "black") +
  theme_bw()


print(plot)



