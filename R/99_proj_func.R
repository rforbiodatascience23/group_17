# Linear model to identify differentially expressed genes.
linear_model_deg <- function(data, test) {
  formula <- as.formula(str_c(test, " ~ value"))
  data |> 
    select(c(test, starts_with("NP"))) |>
    pivot_longer(cols = starts_with("NP"),
                 names_to = "protein",
                 values_to = "value") |>
    group_by(protein) |>
    nest() |> 
    mutate(LM = map(.x = data,
                    .f = ~lm(formula = formula,
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
}

#Volcano plot

volcano_plot <- function(data,type) {
  data <- data |>
    mutate(color_col = case_when(
      adjusted_fdr < 0.05 & estimate > 0 ~ "Significant positive influence", 
      adjusted_fdr < 0.05 & estimate < 0 ~ "Significant negative influence",
      TRUE ~"Not significant"), 
    log2 = -log2(p.value))
  ggplot(data = data, aes(x = estimate, y = log2)) +
    geom_point(aes(color = color_col)) +
    geom_hline(yintercept = 18, linetype = "dashed", color = "black" ) +
    geom_hline(yintercept = (2^-0.8), linetype = "dashed", color = "black" ) +
    scale_color_manual(values = c("Significant positive influence" = "orange", 
                                  "Significant negative influence" = "#69b3a2",
                                  "No significant association" = "grey")) +
    geom_text_repel(data = data |>filter(adjusted_fdr < 0.05),
                    aes(label = protein), hjust = 0, vjust = 0, size = 3, 
                    color = "black") +
    theme_bw() + 
    labs(title = str_c("volcano plot of ", type),
         color = "Significance level ",
              x = "Estimate", 
              y = "-Log2(FDR-adjusted P-value)")
}

#Coorrelation plot

corr_plot <- function(data,test){
  data_select <- data |>
    select(c("id", "er_status", "pr_status", "her2_final_status", test)) |> 
    drop_na()
  
  data_corr <- correlate(data_select, method = "kendall") |>
    focus(er_status, pr_status, her2_final_status)
  
  data_corr <- data_corr |> 
    dplyr::rename(Protein = term) |>   
    pivot_longer(cols = ends_with("status"),
                 names_to = "status", 
                 values_to = "corr") |> 
    mutate(
      status = recode(status, 
                      "er_status" = "Estrogen", 
                      "pr_status" = "Progesteron",
                      "her2_final_status" = "her2"),
      status = factor(status, levels = unique(status)))

  data_corr|>
  ggplot(aes(x = status, y = Protein, fill = corr)) +
    geom_tile(height=0.8, width=0.8) +
    geom_text(aes(label = round(corr,2)), color = "black", size = 3) +
    scale_fill_gradient2(low="blue", mid="white", high="red", limits = c(-1,1)) +
    theme_minimal() +
    coord_equal() +
    theme(plot.title = element_text(hjust = 0.4, size = 12),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1, margin = margin(-3, 0, 0, 0)),
          axis.text.y = element_text(size = 8, margin = margin(0, -3, 0, 0)),
          panel.grid.major = element_blank())
}
