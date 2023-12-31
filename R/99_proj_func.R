##########################################
# Hierarchical clustering and dendrograms 
##########################################

hc_dendro_plot <- function(count_matrix, data_type){
  # create distance matrix (scaled) and perform hierachial clustering
  hc <- count_matrix |>
    stats::dist() |> # calculating the distance for each value
    stats::hclust(method = "ward.D2")  # applying the hierachial clustering algorithm to the data 
  
  # save clustering as dendogram
  ddg <- stats::as.dendrogram(hc)
  
  # extract dendogram order
  ddg_order <- stats::order.dendrogram(ddg)
  
  # extract all dendro information
  ddg_data <- ggdendro::dendro_data(ddg)
  
  # add pam50 information to ddg_data
  ddg_label <- ggdendro::label(ddg_data) |>
    dplyr::rename(id = label) |> 
    dplyr::left_join(clinical_data, by = c("id")) |> 
    dplyr::select(x,y, id, pam50)
  
  # create dendogram plot based on hierachial clustering
  ddg_plot <- ggdendro::ggdendrogram(data = hc, rotate=TRUE) +
    theme(axis.text.y = element_blank(), # removing labels on x axis
          axis.text.x = element_blank(), # removing labels on Y axis
          plot.title = element_text(hjust = 0.5, vjust = 0)) +  
    geom_text(data=ddg_label,
              aes(label=id, x=x, y=-45, colour=as.factor(pam50)), # labels on y axis colored by pam50. 
              size = 1.5) +
    labs(title = str_c("HC based on", data_type, sep = " "))
  return(ddg_plot)
}



########################################################
# Linear model to identify differential expressed genes.
########################################################

linear_model_deg <- function(data, test) {
  formula <- as.formula(str_c(test, " ~ value"))
  data |> 
    select(c(test, starts_with("NP"))) |>
        # Select only the protein data and the column storing information on 
        # whether a person is the HER2 group of not (binary)
    pivot_longer(cols = starts_with("NP"),
                 names_to = "protein",
                 values_to = "value") |>
        # Pivot longer, listing each protein under one another
    group_by(protein) |>
        # Group by protein
    nest() |> 
        # Nest according to the grouping. Data is now nested by proteins. 
    mutate(LM = map(.x = data,
                    .f = ~lm(formula = formula,
                             data = .x))) |>
        # Take one nested dataframe (one protein) at a time and run a linear model 
        # where the protein value is the independent variable, and the HER2-status
        # is the dependent variable.
        # store the result in a column in the mother-dataframe called "LM"
        # The results are now stored in lists nested in the mother dataframe. 
    mutate(LM_tidy = map(.x = LM,
                         .f = ~broom::tidy(x = .x,
                                           conf.int = TRUE,
                                           conf.level = 0.95))) |>
        # Unpack the results from the linear models. Store them in a column 
        # called 'LM_tidy'
    unnest(LM_tidy) |>
       # Unpack the results into the mother dataframe for easier readability
    filter(term == "value") |>
        # Filter the columns to include only "value". This is done because the last
        # line of code gave 2 rows of results for each protein - one for intercept
        # and one for protein. We are only interested in the protein.
    select(-c(data, LM, term)) |>
        # Deselect the column storing the raw protein data to only have the results left
    mutate(adjusted_BH = p.adjust(p.value, n = 12554, method = "BH"), 
               # Benjamini hochberg adjustment for multiple comparisons
           adjusted_fdr = p.adjust(p.value, n = 12554, method = "fdr")) |>
               # False discrovery rate adjustment for multiple comparisons
    arrange(p.value)
               # arrange data according to best adjusted p-values
}

#################
#Volcano plots
#################

volcano_plot <- function(data,subtype) {
  data <- data |>
    mutate(
      color_col = case_when(
        adjusted_fdr < 0.05 & estimate > 0 ~ "Significant positive influence",
        adjusted_fdr < 0.05 & estimate < 0 ~ "Significant negative influence",
        TRUE ~"Not significant"), 
      log2 = -log2(p.value))
      #creating new variable color_col, where the condition case_when checks whether #adjusted fdr is lees than 0.05 and estimate smaller or larger than 0.
      #TRUE-"not signficant is default value of no condition is met."
  
  ggplot(data = data, aes(x = estimate, y = log2)) +
    geom_point(aes(color = color_col)) +
    geom_hline(yintercept = (2^-0.8), linetype = "dashed", color = "black" ) +
    geom_hline(yintercept = 18, linetype = "dashed", color = "black" ) +
    scale_color_manual(values = c("Significant positive influence" = "orange", 
                                  "Significant negative influence" = "#69b3a2",
                                  "No significant association" = "grey")) +
    geom_text_repel(data = data |>filter(adjusted_fdr < 0.05),
                    aes(label = protein), hjust = 0, vjust = 0, size = 3, 
                    color = "black") +
    annotate(geom="text", x=max(data |> pull(estimate)) - 0.09, y=18, label="FDR-adjusted significance level", color="black", size = 3, vjust = 1) +
    theme_bw()+
    labs(title = str_c("Volcano plot of ", subtype),
         color = "Significance level ",
              x = "Estimate", 
              y = "-Log2(FDR-adjusted P-value)")
  
          #Initialize a ggplot with data from the HER2 data frame and aesthetics mappings are #specified.
          #adds points to ploe using color_col varaible for color
          #labs. sets title and labels
          # geom_hline adds horizontal lines to the plot
          #scale_color_manual adds specific color to vlaues dependent of significans levels.
          #geom_text_rebel add labelse to significant points.
          #theme_bw add black and white theme.
          #labs adds title overall and color legends title
          #plot for HER2
}


#####################
#Correlation plots
#####################


Corr_subclasses <- function(data) {
  data|>
    ggplot(aes(x = subclass, y = status, fill = corr)) +
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
          panel.grid.major = element_blank()) + 
    labs(title = test)
}






############################################
# Heatmap of results from diff exp analysis
############################################

# Small function for extracting the neccesary data from the differential expression analysis results and renaming:
subset_rename_data <- function(data, subtype) {
  cols = c("estimate", "adjustedfdr")
  data <- data |>
    dplyr::rename("adjustedfdr" = "adjusted_fdr") |>
    select(c("protein", "estimate", "adjustedfdr")) |>
    dplyr::rename_with(.fn = ~paste0(.,"_",subtype), any_of(cols))
}

# Heatmaps:

heat_map_results <- function(data, fdr_col, title){
  data_long <- data |>
    arrange({{fdr_col}}) |>
    slice(1:20)|>
    # Select only the top 20 most significant proteins for subtype of interest
    mutate(protein = factor(protein), 
           protein = fct_reorder(protein, {{fdr_col}})) |>
    # Make sure protein is a factor varaible + make the levels the order in 
    # which each prtein appears. This is for vizualization: if levels are not specified,
    # the heatmap plots the proteins in alphabetical order. We want to plot
    # them in order of significance, and thus, determining the levels overrules the 
    # alphabetical order. 
    pivot_longer(cols = -protein,
                 names_to = c(".value", "pam50"),
                 names_sep = "_") |>
    drop_na(pam50) |> # removing columns that contain na in the pam50 column.
    # Pivot longer - because this is how geom_heatmap() reads the data
    mutate(
      pam50 = case_when(
        pam50 == "HER2" ~ "HER2 enriched", 
        pam50 == "lumA" ~ "Luminal A", 
        pam50 == "lumB" ~ "Luminal B", 
        pam50 == "bl" ~ "Basal like"),
      # Mutating the pam50 column to whole sentences. 
      # Makes the heatmap more readable and prettier
      significance = case_when(
        adjustedfdr < 0.001 ~ "***", 
        adjustedfdr < 0.01 ~ "**",
        adjustedfdr < 0.05 ~ "*",
        TRUE ~ ""
      )
      # Add a column of asterisks showing significance. Also to be used in vizualisarion
    ) 

  
  heatmap <- ggplot(data_long, aes(x = pam50 , y = factor(protein), fill = estimate))+
    # Make plot with pam50 subtypes on x-axis and proteins on y-axis, 
    # Color for the estimate (slope)
    geom_tile() +
    # Add tiles
    scale_fill_gradient2(midpoint = 0, low = "blue", high = "red", mid = "white") +
    # Specify colouring: midpoint (neutral) should be white, and upregulated 
    # protein rea, and downregulated blue
    labs(
      title = title, 
      y = "", 
      x = "", 
      fill = "Estimate"
    ) +
    # Enter labels 
    geom_text(aes(x =pam50 , y = protein, label = significance), size = 4, vjust = 0.8) +
    # Add asterisks describing the significane. This uses the 'significance' column
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          axis.text.y = element_text(size = 5), 
          panel.background = element_blank(),
          plot.title = element_text(size=12))
  # Small adjustements of size, angles and positioning of the text elements in the
  # plot, as well as removal of the grey background. 
  
  return(heatmap)
  #display
}


# In summary, this function takes the dataframe that we use, sort according to significance of the subtype we want to use, and finally subset to only the top 20 most significant protein
# This data-frame is then wrangled - making it into long format + adding columns used for plotting. 
# Finally, the heatmap in itself is plotted.


