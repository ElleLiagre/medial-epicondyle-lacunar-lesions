############################################################################
# Title: Microtrauma at the Humeral Medial Epicondyle - Reliability Study
# 
# Author: Elle B. K. Liagre
# Email: elle.liagre@u-bordeaux.fr
# ORCID: https://orcid.org/0000-0002-8993-3266
# Date: 2025-12-18
# 
# Description: 
#      Reliability analysis for the protocol extracting lacunar lesions from the 
#      humeral medial epicondyle as 3D meshes. Compares surface area measurements
#      between two independent extraction sessions (Set A and Set B) using Lin's 
#      Concordance Correlation Coefficient, with separate analyses for all surfaces
#      and non-marginal surfaces only.
#
# Requirements:
#    - R (>= 4.3.2)
#    - Rvcg (>= 0.25)
#    - dplyr (>= 1.1.4) 
#    - tidyr (>= 1.3.1)
#    - stringr (>= 1.5.1)
#    - ggplot2 (>= 3.5.2)
#    - openxlsx (>= 4.2.5.2)
#    - readxl (>= 1.4.5)
#    - DescTools (>= 0.99.60)
#
# Input data:
#    - 3D mesh files in PLY format from two extraction sessions (Set A and Set B)
#    - Manual mapping file (matching.xlsx) for matching observations between sessions
#
# License: This script is licensed under the GNU General Public License v3.0.
# See https://www.gnu.org/licenses/gpl-3.0.html for more details.
#
###########################################################################

# Dependencies ------------------------------------------------------------

library(Rvcg)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(openxlsx)
library(readxl)
library(DescTools)

# Functions ---------------------------------------------------------------

# Summarise area
summarise_meshes <- function(mesh_list) {
  data.frame(
    mesh_name = names(mesh_list),
    area = sapply(mesh_list, vcgArea),
    stringsAsFactors = FALSE
  )
}

# Import data -------------------------------------------------------------

## Configuration
wd <- getwd()

seta_path <- file.path(wd, "Set_A")
setb_path <- file.path(wd, "Set_B")
matching_file <- file.path(wd, "matching.xlsx")


## Importing of data
# Load Set A meshes
seta_list <- list.files(path = seta_path, pattern = "\\.ply$", full.names = TRUE)

seta <- lapply(seta_list, function(file) {
  vcgPlyRead(file)
})

names(seta) <- tools::file_path_sans_ext(basename(seta_list))
names_seta <- names(seta)


# Load Set B meshes
setb_list <- list.files(path = setb_path, pattern = "\\.ply$", full.names = TRUE)

setb <- lapply(setb_list, function(file) {
  vcgPlyRead(file)
})

names(setb) <- tools::file_path_sans_ext(basename(setb_list))
names_setb <- names(setb)


# Load manual matching between sets
manual_map <- read_excel(matching_file, sheet = "Feuil1")


# Quantification of extracted surfaces ------------------------------------

## Extract mesh properties
seta_df <- summarise_meshes(seta)
setb_df <- summarise_meshes(setb)


## Parse specimen information from filenames
seta_df <- seta_df %>%
  mutate(
    site = str_match(mesh_name, "^AO_([A-Za-z]+)")[,2],
    specimen = str_match(mesh_name, paste0("^AO_[A-Za-z]+_(.+?)_remeshed"))[,2],
    comp = str_remove(str_extract(mesh_name, "comp[0-9_]+"), "comp")
  )


setb_df <- setb_df %>%
  mutate(
    site = str_match(mesh_name, "^AO_([A-Za-z]+)")[,2],
    specimen = str_match(mesh_name, paste0("^AO_[A-Za-z]+_(.+?)_remeshed"))[,2],
    comp = str_remove(str_extract(mesh_name, "comp[0-9_]+"), "comp")
  )


# Match observations ------------------------------------------------------

## Verify if no mistakes were made in the manual mapping
# Check Set A
map_combinations <- manual_map %>%
  dplyr::select(site, specimen, comp = comp_A) %>%
  dplyr::filter(!is.na(comp) & comp != "NA") %>%  
  dplyr::distinct()

seta_combinations <- seta_df %>%
  dplyr::select(site, specimen, comp) %>%
  dplyr::distinct()

# In manual_map but not in seta_df
missing_in_seta <- anti_join(map_combinations, seta_combinations, by = c("site", "specimen", "comp"))

# In seta_df but not in manual_map
extra_in_seta <- anti_join(seta_combinations, map_combinations, by = c("site", "specimen", "comp"))

# Report
cat("Set A verification:\n")
cat("  Missing in Set A:", nrow(missing_in_seta), "\n")
cat("  Extra in Set A:", nrow(extra_in_seta), "\n")


# Check Set B
map_combinations <- manual_map %>%
  dplyr::select(site, specimen, comp = comp_B) %>%
  dplyr::filter(!is.na(comp) & comp != "NA") %>%  
  dplyr::distinct()

setb_combinations <- setb_df %>%
  dplyr::select(site, specimen, comp) %>%
  dplyr::distinct()

# In manual_map but not in seta_df
missing_in_setb <- anti_join(map_combinations, setb_combinations, by = c("site", "specimen", "comp"))

# In seta_df but not in manual_map
extra_in_setb <- anti_join(setb_combinations, map_combinations, by = c("site", "specimen", "comp"))

# Report
cat("\nSet B verification:\n")
cat("  Missing in Set B:", nrow(missing_in_setb), "\n")
cat("  Extra in Set B:", nrow(extra_in_setb), "\n")


# Combine matched data ----------------------------------------------------

# Combine data in dataframe
matched_df <- manual_map %>%
  left_join(
    seta_df %>% dplyr::select(site, specimen, comp, area),
    by = c("site", "specimen", "comp_A" = "comp")
  ) %>%
  rename(area_A = area) %>%
  left_join(
    setb_df %>% dplyr::select(site, specimen, comp, area),
    by = c("site", "specimen", "comp_B" = "comp")
  ) %>%
  rename(area_B = area)

# Replace all NA with 0
matched_df <- matched_df %>%
  mutate(
    area_A = ifelse(is.na(area_A), 0, area_A),
    area_B = ifelse(is.na(area_B), 0, area_B)
  )

# Export matched data
wb <- createWorkbook()

addWorksheet(wb, "Matched")
writeData(wb, "Matched", matched_df)

saveWorkbook(wb, file = file.path(wd,"matched_data.xlsx"), overwrite = TRUE)


# Descriptive statistics --------------------------------------------------

## Prepare data for stats

stats <- matched_df

# Add sample name (dev / val)
stats$Sample <- ifelse(grepl("^(Bur|Mau|Scl)", stats$site), 
                       "Validation", "Development")

# Map site codes to full names
site_code <- substr(stats$site, 1, 3)
site_map <- c(Bur = "Burnot", 
              Scl = "Sclaigneaux", 
              Mau = "Maurenne", 
              Isl = "Isle Adam", 
              Pet = "Petit Morin")

stats$site <- site_map[site_code]

# Re-order columns and calculate differences
stats <- stats %>%
  dplyr::select(Sample, site, specimen, comp_A, comp_B, area_A, area_B, everything())

stats$diff_area <- abs(stats$area_A - stats$area_B)


## Calculate descriptive statistics

desc_stats <- stats %>%
  summarise(
    across(
      c(diff_area),
      list(
        mean       = ~mean(.x, na.rm = TRUE),
        sd         = ~sd(.x, na.rm = TRUE),
        min        = ~min(.x, na.rm = TRUE),
        q1         = ~quantile(.x, 0.25, na.rm = TRUE),
        median     = ~median(.x, na.rm = TRUE),
        q3         = ~quantile(.x, 0.75, na.rm = TRUE),
        max        = ~max(.x, na.rm = TRUE),
        n_missing  = ~sum(is.na(.x)),
        n          = ~sum(!is.na(.x))
      ),
      .names = "{.col}.{.fn}"
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", "Statistic"),
    names_sep = "\\."
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = value
  )


## Visualisation of descriptive statistics

diff_long <- stats %>%
  dplyr::select(diff_area) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Difference")

# Boxplot
box_area <- ggplot(diff_long %>% filter(Variable == "diff_area"), aes(x = "", y = Difference)) +
  geom_boxplot() +
  labs(
    title = "Absolute Differences for Area",
    y = "Absolute Difference (mm²)",
    x = ""
  ) +
  theme_minimal(base_size = 14)

# Histogram
hist_area <- ggplot(diff_long %>% filter(Variable == "diff_area"), aes(x = Difference)) +
  geom_histogram() +
  labs(
    title = "Histogram of Absolute Differences for Area",
    x = "Absolute Difference (mm²)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

# Scatter plot
scatter_area <- ggplot(stats, aes(x = area_A, y = area_B)) +
  geom_point() +
  labs(
    title = "Scatter Plot of Surface Areas",
    x = "Surfaces of Set A (mm²)",
    y = "Surfaces of Set B (mm²)"
  ) +
  theme_minimal(base_size = 14)


## Export results

wb <- createWorkbook()

addWorksheet(wb, "Stats")
writeData(wb, "Stats", stats %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

addWorksheet(wb, "Descr_Stats")
writeData(wb, "Descr_Stats", desc_stats %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

ggsave(file.path(wd, "boxplot_area.png"), box_area, width = 6, height = 4, dpi = 600)
ggsave(file.path(wd, "histogram_area.png"), hist_area, width = 6, height = 4, dpi = 600)
ggsave(file.path(wd, "scatter_area.png"), scatter_area, width = 6, height = 4, dpi = 600)


# Lin's CCC ---------------------------------------------------------------

## Calculation
area_ccc <- CCC(stats$area_A, stats$area_B)


ccc_df <- data.frame(
  variable = "Area",
  rho_c = area_ccc$rho.c$est,
  lwr_ci = area_ccc$rho.c$lwr.ci,
  upr_ci = area_ccc$rho.c$upr.ci
)


## Export results

addWorksheet(wb, "LinCCC")
writeData(wb, "LinCCC", ccc_df %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

saveWorkbook(wb, file = file.path(wd,"stats_results.xlsx"), overwrite = TRUE)


# Statistical tests (no margin) -------------------------------------------

stats_no_margin <- stats %>%
  filter(cropping_margin != 1)


## Descriptive statistics

desc_stats_no_margin <- stats_no_margin %>%
  summarise(
    across(
      c(diff_area),
      list(
        mean       = ~mean(.x, na.rm = TRUE),
        sd         = ~sd(.x, na.rm = TRUE),
        min        = ~min(.x, na.rm = TRUE),
        q1         = ~quantile(.x, 0.25, na.rm = TRUE),
        median     = ~median(.x, na.rm = TRUE),
        q3         = ~quantile(.x, 0.75, na.rm = TRUE),
        max        = ~max(.x, na.rm = TRUE),
        n_missing  = ~sum(is.na(.x)),
        n          = ~sum(!is.na(.x))
      ),
      .names = "{.col}.{.fn}"
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", "Statistic"),
    names_sep = "\\."
  ) %>%
  pivot_wider(
    names_from = Statistic,
    values_from = value
  )


## Visualisation of descriptive statistics

# Gather the two difference columns into long format for easier plotting
diff_long <- stats_no_margin %>%
  dplyr::select(diff_area) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Difference")

# Boxplot
box_area <- ggplot(diff_long %>% filter(Variable == "diff_area"), aes(x = "", y = Difference)) +
  geom_boxplot() +
  labs(
    title = "Absolute Differences for Area",
    y = "Absolute Difference (mm²)",
    x = ""
  ) +
  theme_minimal(base_size = 14)


# Histogram
hist_area <- ggplot(diff_long %>% filter(Variable == "diff_area"), aes(x = Difference)) +
  geom_histogram() +
  labs(
    title = "Histogram of Absolute Differences for Area",
    x = "Absolute Difference (mm²)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)


# Scatter plot
scatter_area <- ggplot(stats_no_margin, aes(x = area_A, y = area_B)) +
  geom_point() +
  labs(
    title = "Scatter Plot of Surface Areas",
    x = "Surfaces of Set A (mm²)",
    y = "Surfaces of Set B (mm²)"
  ) +
  theme_minimal(base_size = 14)


## Export results

wb <- createWorkbook()

addWorksheet(wb, "Stats_no_margin")
writeData(wb, "Stats_no_margin", stats_no_margin %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

addWorksheet(wb, "Descr_Stats_no_margin")
writeData(wb, "Descr_Stats_no_margin", desc_stats_no_margin %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

ggsave(file.path(wd, "boxplot_area_no_margin.png"), box_area, width = 6, height = 4, dpi = 600)
ggsave(file.path(wd, "histogram_area_no_margin.png"), hist_area, width = 6, height = 4, dpi = 600)
ggsave(file.path(wd, "scatter_area_no_margin.png"), scatter_area, width = 6, height = 4, dpi = 600)



### Lin's CCC
area_ccc <- CCC(stats_no_margin$area_A, stats_no_margin$area_B)


ccc_df <- data.frame(
  variable = "Area",
  rho_c = area_ccc$rho.c$est,
  lwr_ci = area_ccc$rho.c$lwr.ci,
  upr_ci = area_ccc$rho.c$upr.ci
)


## Export results

addWorksheet(wb, "LinCCC")
writeData(wb, "LinCCC", ccc_df %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

saveWorkbook(wb, file = file.path(wd,"stats_no_margin_results.xlsx"), overwrite = TRUE)

