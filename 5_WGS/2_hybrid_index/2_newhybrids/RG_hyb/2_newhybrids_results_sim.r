#Analyzing NewHybrid output

# --- 0. Load necessary libraries ---
# Ensure these libraries are installed with install.packages() if you haven't already.
.libPaths(c(.libPaths(),"/dss/dsshome1/lxc0B/ra52qed/R/x86_64-pc-linux-gnu-library/4.2"))
library(dplyr)
library(tidyr)
library(stringr)
library(adegenet)
library(ggplot2)
library(hybriddetective)
library(genepopedit)
library(parallelnewhybrid)

# ---- 1. Read in simulation results -----------
# Set your path to the simulation results
sim_path <- "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RG_hyb/1_simulations/NH.Results/"

# Find all PofZ files in all subfolders
sim_files <- list.files(path = sim_path, 
                        pattern = "_PofZ.txt$", 
                        full.names = TRUE, 
                        recursive = TRUE)

all_sims_list <- NULL

for (f in sim_files) {
  # NewHybrids output has no header, or generic headers
  temp_dat <- read.table(f, header = TRUE) 
  
  # Rename columns to standard classes using dplyr
  temp_dat <- temp_dat %>%
    dplyr::rename(P0 = 3, P1 = 4, F1 = 5, F2 = 6, BC1 = 7, BC2 = 8)
  
  # Extract metadata from the folder path
  path_parts <- unlist(strsplit(f, "/"))
  rep_name <- path_parts[length(path_parts)-1]
  
  # Determine "Known" Ancestry
  # We assign known types to rows
  temp_dat <- temp_dat %>%
    mutate(Replicate = rep_name,
           Known_Ancestry = case_when(
             row_number() %in% 1:12  ~ "Pure_Pop1",
             row_number() %in% 13:24 ~ "Pure_Pop2",
             row_number() %in% 25:39 ~ "F1",
             row_number() %in% 40:54 ~ "F2",
             row_number() %in% 55:69 ~ "BC1",
             row_number() %in% 70:84 ~ "BC2",
             TRUE ~ "Other"
           ))
  
  all_sims_list[[f]] <- temp_dat
}

# Merge all simulations into one dataframe for Power Analysis
all_sim_results <- bind_rows(all_sims_list)


# ----- 2. Assign categories based on your chosen threshold-------
# combine data from all 9 simulation runs
sim_assigned <- all_sim_results %>%
  mutate(Assigned_Category = case_when(
    P0 >= 0.995 ~ "Pure_Pop1",
    P1 >= 0.995 ~ "Pure_Pop2",
    F1 >= 0.995 ~ "F1",
    F2 >= 0.995 ~ "F2",
    BC1 >= 0.995 ~ "BC1",
    BC2 >= 0.995 ~ "BC2",
    TRUE ~ "Indeterminate"
  ))

# 2. Calculate Efficiency and Accuracy for each class
performance_metrics <- sim_assigned %>%
  group_by(Known_Ancestry) %>%
  summarize(
    # Efficiency: (True Positives) / (Total Known)
    Efficiency = sum(Assigned_Category == Known_Ancestry) / n(),
    
    # Accuracy: (True Positives) / (Total Assigned to that category)
    # This requires looking at the whole dataset for each class
    .groups = 'drop'
  )

# 3. For Accuracy, we need a different grouping
accuracy_metrics <- sim_assigned %>%
  filter(Assigned_Category != "Indeterminate") %>%
  group_by(Assigned_Category) %>%
  summarize(
    Accuracy = sum(Assigned_Category == Known_Ancestry) / n(),
    .groups = 'drop'
  )

# ---- 3. Run this in a loop with many cut-offs ----
thresholds <- seq(0.8, 0.995, by = 0.01)
results_list <- list()

for(t in thresholds){
  
  # 1. Assign categories based on current threshold 't'
  sim_temp <- all_sim_results %>%
    mutate(Assigned = case_when(
      P0 >= t ~ "Pure_Pop1",
      P1 >= t ~ "Pure_Pop2",
      F1 >= t ~ "F1",
      F2 >= t ~ "F2",
      BC1 >= t ~ "BC1",
      BC2 >= t ~ "BC2",
      TRUE ~ "Indeterminate"
    ))
  
  # 2. Calculate Efficiency per class (True Positives / Total Known)
  eff <- sim_temp %>%
    group_by(Known_Ancestry) %>%
    summarise(Efficiency = sum(Assigned == Known_Ancestry) / n(), .groups = 'drop')
  
  # 3. Calculate Accuracy per class (True Positives / Total Assigned)
  acc <- sim_temp %>%
    filter(Assigned != "Indeterminate") %>%
    group_by(Assigned) %>%
    summarise(Accuracy = sum(Assigned == Known_Ancestry) / n(), .groups = 'drop') %>%
    dplyr::rename(Known_Ancestry = Assigned) # Rename for merging
  
  # 4. Merge and store
  perf <- full_join(eff, acc, by = "Known_Ancestry") %>%
    mutate(Threshold = t)
  
  results_list[[as.character(t)]] <- perf
}

# Combine all results into one dataframe
df_performance <- bind_rows(results_list) %>%
  mutate(Overall_Performance = Efficiency * Accuracy)


#Plot overall performance
ggplot(df_performance, aes(x = Known_Ancestry, y = Overall_Performance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Hybrid Class",
    y = "Overall Performance",
    title = "NewHybrids Performance Across Simulation Iterations"
  )


# To see which threshold gives the best "Overall Performance" (Eff * Acc)
best_thresholds <- df_performance %>%
  group_by(Known_Ancestry) %>%
  filter(Overall_Performance == max(Overall_Performance, na.rm = TRUE))

print(best_thresholds)

##Assignment of F1 and Pure parentals above 0.6
##Advanced hybrids low reliability

#### ---- Empirical Results --------

emp_poz <- read.table(
  "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RG_hyb/2_empirical/NH.Results/RG_70_NH.txt_Results/RG_70_NH.txt_PofZ.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

ind_poz <- read.table(
  "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RG_hyb/2_empirical/NH.Results/RG_70_NH.txt_Results/RG_70_NH_individuals.txt",
  header = F,
  stringsAsFactors = FALSE
)

emp_poz$Individual <- ind_poz$V1
emp_poz <- emp_poz %>%
  dplyr::rename(P0 = 3, P1 = 4, F1 = 5, F2 = 6, BC1 = 7, BC2 = 8)

emp_assigned <- emp_poz %>%
  mutate(
    Assigned_Class = case_when(
      P0  >= 0.90 ~ "G",
      P1  >= 0.90 ~ "R",
      F1  >= 0.90 ~ "F1",
      F2  >= 0.90 ~ "F2",
      BC1 >= 0.90 ~ "BC_G",
      BC2 >= 0.90 ~ "BC_R",
      TRUE        ~ "Indeterminate"
    )
  )


write.table(emp_assigned, file= "/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/5_WGS/2_hybrid_index/2_newhybrids/RG_hyb/2_output_classes.txt", quote=F, row.names=F, sep="\t")
