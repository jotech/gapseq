# Weights distribution plot
rxn.weights.dist.plot <- function(weights_plot_dt, dummy.weight, output.dir) {
    num.mod <- length(unique(weights_dt$model_id))
    name.rxn <- unique(weights_plot_dt$seed)
    num.rxn <- length(name.rxn) 
    weights_plot_dt <- weights_plot_dt[,c("seed", "weight")]

    # Add a new column with a unique ID composed of seed + progressive number
    weights_plot_dt <- weights_plot_dt %>%
    group_by(seed) %>%
    mutate(unique_id = paste0(seed, "_", row_number())) %>%
    ungroup() %>%
    select(-seed) 

    # Create a data.frame with all possible seeds-weight
    seed = rep(name.rxn, each = num.mod)
    seq.tmp = rep(seq(1,num.mod),num.rxn)
    unique_id <- paste(seed, seq.tmp, sep = "_")
    all_seeds <- data.frame(
    seed = seed,
    unique_id = unique_id)

    # Merge the weight dt and the complete empty df to obtain a rxn weight-score in each MAG
    result_df <- full_join(weights_plot_dt, all_seeds, by = "unique_id") # Merge the existing data frame with all possible seeds
    result_df$weight[is.na(result_df$weight)] <- dummy.weight  # Replace missing weight values with 100
    result_df <- result_df %>% 
    arrange(seed, weight) %>%
    mutate(x.axis = rep(seq(1,num.mod),num.rxn))

    # weight_score_dist --- PLOT
    p <- ggplot(result_df, aes(x = x.axis, y = weight, color = factor(seed))) +  
        geom_line() +
        theme_minimal() +
        labs(
        x = "MAG",
        y = "Weight") +
        guides(color = "none")  # Turn off the legend for the "seed" variable
    ggsave(file.path(output.dir,"weight_score_dist.pdf"), p, width = 8, height = 6, units = "in", dpi = 300)

    # cumulative_weight_score_dist --- PLOT
    p_df <- result_df %>%
    group_by(x.axis) %>%
    summarise(sum_weight = sum(weight)/num.rxn)
    
    p <- ggplot(p_df, aes(x = x.axis, y = sum_weight)) +  
        geom_bar(stat = "identity", fill = "blue") + 
        theme_minimal() +
        labs(
        x = "MAG",
        y = "Weight") 
    ggsave(file.path(output.dir,"cumulative_weight_score_dist.pdf"), p, width = 8, height = 6, units = "in", dpi = 300)
}

# RXN frequency and U-shape plot
rxn.freq.plot <- function(rxn2mod_dt, output.dir) {
    num.mod <- dim(rxn2mod_dt)[2]-1 
    # How many times a RXN have been detected
    rxn_freq_df <- data.frame(rxn = rxn2mod_dt$rxn, frequency = rowSums(rxn2mod_dt[,-"rxn"])) 
    rxn_freq_df <- rxn_freq_df[order(rxn_freq_df$frequency, rxn_freq_df$rxn), ] # sort by frequency and alphabetically
    # How many RXN have a certain frequency
    num_detected_fred_df <- data.frame(rxn_frequency = names(table(rxn_freq_df$frequency)), num_rxn = as.vector(table(rxn_freq_df$frequency))) # some frequencies will never be detected so they will be missing
    tmp_df <- data.frame(rxn_frequency = seq(1, num.mod), num_rxn = rep(0, num.mod)) # empty df to include missing RXN frequencies
    tmp_df <- merge(tmp_df, num_detected_fred_df, by = "rxn_frequency", all = TRUE)
    tmp_df$num_rxn <- ifelse(is.na(tmp_df$num_rxn.y), tmp_df$num_rxn.x, tmp_df$num_rxn.y) # if RXN frequency was missing add 0
    num_fred_df <- tmp_df[, c("rxn_frequency", "num_rxn")]
    num_fred_df$rxn_frequency <- factor(num_fred_df$rxn_frequency, levels = num_fred_df$rxn_frequency)

    # rxn_frequency_plot --- PLOT
    x_labels <- list()
    for (i in seq(1, length(rxn_freq_df$rxn), by = 30)) { # plot only some ticks in the x-axis
    x_labels <- c(x_labels, as.character(rxn_freq_df$rxn[i]))
    }
    rxn_freq_df$rxn <- factor(rxn_freq_df$rxn,                                    
                    levels = rxn_freq_df$rxn)
    rxn_frequency_plot <- ggplot(data=rxn_freq_df, aes(x=rxn, y=frequency)) + 
        geom_bar(stat="identity", fill="grey",  width=1)+
        scale_x_discrete(breaks = unlist(x_labels))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(file.path(output.dir,"rxn_frequency_plot.pdf"), rxn_frequency_plot, width = 14, height = 6, units = "in", dpi = 300)
    
    # u_shape_plot --- PLOT
    x_labels <- list() # plot only some ticks in the x-axis
    for (i in seq(1, length(num_fred_df$rxn_frequency), by = 20)) {
    x_labels <- c(x_labels, as.character(num_fred_df$rxn_frequency[i]))
    }
    u_shape_plot <- ggplot(data=num_fred_df, aes(x=rxn_frequency, y=num_rxn)) +    
        geom_bar(position=position_dodge(0.9), stat="identity", fill="grey",  width=0.9)+
        scale_x_discrete(breaks = unlist(x_labels))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"))
    ggsave(file.path(output.dir,"u_shape_plot.pdf"), u_shape_plot, width = 8, height = 6, units = "in", dpi = 300)
}

### Rarefaction courve plot
rarefaction_mean_std <- function(dt, curve_type) {
  means <- t(dt[, lapply(.SD, mean, na.rm = TRUE)])
  colnames(means) <- "MEAN"
  std_devs <- t(dt[, lapply(.SD, sd, na.rm = TRUE)])
  colnames(std_devs) <- "STD"
  single_rarefaction_dt <- data.table(means, std_devs) #, by.x = "MEAN", by.y = "STD", all = TRUE)
  single_rarefaction_dt[, curve_type := rep(curve_type, dim(single_rarefaction_dt)[1])]
  single_rarefaction_dt[, num_models := 1:dim(single_rarefaction_dt)[1]]

  return(single_rarefaction_dt)
}

rarefaction_curve_data <- function(rarefaction_dt, num_itreation, core_th) {
  cat("Iterating models to generate rarefaction data... \n")

  num_mods <- dim(rarefaction_dt)[2]
  cum_num_rxn_df <- matrix(NA, nrow = num_itreation, ncol = num_mods)
  cum_num_rxn_CORE_df <- matrix(NA, nrow = num_itreation, ncol = num_mods)
  cum_num_rxn_ACC_df <- matrix(NA, nrow = num_itreation, ncol = num_mods)

  for (iter in 1:num_itreation){
    # cat(paste(iter, "\n"))
    # randomly shuffle models names
    randomly_shuffle_models <- sample(colnames(rarefaction_dt))

    for (idx in seq_along(randomly_shuffle_models)){
      columns_to_extract <- randomly_shuffle_models[1:idx]
      extracted_columns <- rarefaction_dt[, ..columns_to_extract]
      row_sums <- rowSums(extracted_columns)
      # determine the value for the ALL
      num_rxn_all <- sum(row_sums >= 1)
      # determine the value for the CORE and ACCESSORY 
      num_rxn_core <- sum(row_sums >= round(idx*core_th))
      num_rxn_accessory <- sum(row_sums < round(idx*core_th))

      cum_num_rxn_df[[iter, idx]] <- num_rxn_all
      cum_num_rxn_CORE_df[[iter, idx]] <- num_rxn_core
      cum_num_rxn_ACC_df[[iter, idx]] <- num_rxn_accessory
    }
    cum_num_rxn_ACC_df[[iter, 1]] <- 0
  }
  
  # Set column names to the generated numbers
  cum_num_rxn_dt <- data.table(cum_num_rxn_df)
  cum_num_rxn_CORE_dt <- data.table(cum_num_rxn_CORE_df)
  cum_num_rxn_ACC_dt <- data.table(cum_num_rxn_ACC_df)

  # set colnames
  setnames(cum_num_rxn_dt, colnames(cum_num_rxn_dt), as.character(1:length(randomly_shuffle_models)))
  setnames(cum_num_rxn_CORE_dt, colnames(cum_num_rxn_CORE_dt), as.character(1:length(randomly_shuffle_models)))
  setnames(cum_num_rxn_ACC_dt, colnames(cum_num_rxn_ACC_dt), as.character(1:length(randomly_shuffle_models)))
  
  # generate univocous data.table
  rarefaction_dt <- rarefaction_mean_std(cum_num_rxn_dt, "All")
  rarefaction_dt <- rbind(rarefaction_dt, rarefaction_mean_std(cum_num_rxn_CORE_dt, "Core"))
  rarefaction_dt <- rbind(rarefaction_dt, rarefaction_mean_std(cum_num_rxn_ACC_dt, "Shell"))

  return(rarefaction_dt)
}

# Given a dataset with presence absence of reaction in each genome compute the rarefaction courve
#   - X_axis (number of genome considered)
#   - Y_axis (number of unique reactions in the pan-Draft)
#   - STD of the number of unique reactions in the model
#   - Category of the rxn: "Core" or "all"
rarefaction.curve.plot <- function(rxn2mod_dt, output.dir, core_th, num_iter, script.dir) {
    rxn2mod_rarefaction_dt <- rxn2mod_dt[, .SD, .SDcols = setdiff(names(rxn2mod_dt), c("rxn"))]
    rownames(rxn2mod_rarefaction_dt) <- rxn2mod_dt$rxn
    rxn2mod_rarefaction_dt[, (names(rxn2mod_rarefaction_dt)) := lapply(.SD, as.logical), .SDcols = names(rxn2mod_rarefaction_dt)] # Iterate through columns and convert integers to boolean
    rarefaction_dt <- rarefaction_curve_data(rxn2mod_rarefaction_dt, num_iter, core_th)

    # rarefaction_curve_plot --- PLOT
    rarefaction_p <- ggplot(rarefaction_dt, aes(x=num_models, y=rarefaction_dt$MEAN, group=curve_type, color=curve_type)) + 
        geom_line() +
        geom_point()+
        geom_errorbar(aes(ymin=MEAN-STD, ymax=MEAN+STD), width=.2,
                        position=position_dodge(0.05))+
        labs(x="Number of Draft", y = "Mean number of RXN")+
        theme(   
            #plot.background = element_rect(fill = "white"),
            #panel.background = element_rect(fill = "white"),
            aspect.ratio = 1/1.3,  # Adjust the aspect ratio to make x-axis longer
            plot.margin = margin(10, 30, 10, 10)  # Adjust margins for the plot) +
            )+
            scale_color_manual(values=c('#2F8CCA','#CA3D2F','#3DB437'))
    ggsave(file.path(output.dir,"rarefaction_curve_plot.pdf"), rarefaction_p, width = 14, height = 6, units = "in", dpi = 300)
}


