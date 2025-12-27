############################################################################
# Title: Microtrauma at the Humeral Medial Epicondyle - Lesion Analysis
# 
# Author: Elle B. K. Liagre
# Email: elle.liagre@u-bordeaux.fr
# ORCID: https://orcid.org/0000-0002-8993-3266
# Date: 2025-12-18
# 
# Description: 
#   This script analyzes 3D surface meshes of lesions, extracted from cropped
#   medial epicondyles. It computes geometric metrics (depth, volume, surface 
#   area), shape descriptors (sphericity, elongation, flatness), and surface
#   characteristics (solidity, roughness) from PLY mesh files. The analysis
#   includes alignment of meshes, curvature mapping, and statistical summaries.
#
# Requirements:
#    - R (>= 4.3.2)
#    - dplyr (>= 1.1.4)
#    - doolkit (>= 1.42.2)
#    - FNN (>= 1.1.4.1)
#    - GGally (>= 2.2.1)
#    - Hmisc (>= 5.2.4)
#    - openxlsx (>= 4.2.5.2)
#    - reshape2 (>= 1.4.4)
#    - Rvcg (>= 0.25)
#    - tibble (>= 3.2.1)
#    - vangogh (>= 0.1.2.)
#
# Input data:
#     - Extracted lesion surface components as PLY files
#     - Corresponding cropped humeri as PLY files
#     - Capped surface meshes as PLY files
#     - Volume calculations of capped surfaces from PyMeshLab (CSV file)
#     - Volume calculations of convex hulls from PyMeshLab (CSV file)
#
# License: This script is licensed under the GNU General Public License v3.0.
# See https://www.gnu.org/licenses/gpl-3.0.html for more details.
#
###########################################################################


# Clean workspace

rm(list=ls())


# Dependencies ------------------------------------------------------------


library(dplyr)
library(doolkit)
library(FNN)
library(GGally)
library(Hmisc)
library(openxlsx)
library(reshape2)
library(Rvcg)
library(tibble)
library(vangogh)


# Functions ---------------------------------------------------------------

### Mesh manipulation functions

## Align mesh to boundary plane
align_to_boundary_plane <- function(mesh, export_path = NULL, 
                                    format = c("ply", "obj", "stl")) {
  format <- match.arg(format)
  
  # Extract vertices
  vertices <- t(mesh$vb[1:3, ])
  
  # Get boundary vertices
  edges <- Rvcg::vcgGetEdge(mesh)
  boundary_idx <- unique(c(edges$vert1[edges$border == 1], 
                           edges$vert2[edges$border == 1]))
  if (length(boundary_idx) == 0) boundary_idx <- 1:nrow(vertices)
  boundary_vertices <- vertices[boundary_idx, ]
  
  # Fit PCA to boundary vertices to determine principal plane
  pca <- prcomp(boundary_vertices, center = TRUE)
  normal <- pca$rotation[, 3]
  centroid <- colMeans(boundary_vertices)
  
  # Center vertices at origin
  vertices_centered <- sweep(vertices, 2, centroid, "-")
  
  # Compute rotation matrix to align normal to Z-axis
  target <- c(0, 0, 1)
  v <- pracma::cross(normal, target)
  s <- sqrt(sum(v^2))
  
  if (s < 1e-6) {

    # Normal already aligned with Z-axis
    vertices_aligned <- vertices_centered

  } else {
    
    # Rodrigues' rotation formula
    c_val <- sum(normal * target)
    vx <- matrix(c(
      0, -v[3], v[2],
      v[3], 0, -v[1],
      -v[2], v[1], 0
    ), 3, 3, byrow = TRUE)
    R <- diag(3) + vx + (vx %*% vx) * ((1 - c_val) / s^2)
    vertices_aligned <- t(R %*% t(vertices_centered))
  }
  
  # Check orientation: if majority of vertices are above Z=0, flip
  z_vals <- vertices_aligned[, 3]
  pos_ratio <- mean(z_vals > 0)
  
  if (pos_ratio > 0.5) {
    vertices_aligned[, 3] <- -vertices_aligned[, 3]
    message("Mesh flipped along Z-axis (normal inverted).")
  }
  
  # Update mesh vertices
  mesh$vb[1:3, ] <- t(vertices_aligned)
  mesh$it <- mesh$it[c(1,3,2), ]
  mesh <- Rvcg::vcgUpdateNormals(mesh, silent = TRUE)
  
  
  # Optional export
  if (!is.null(export_path)) {
    if (format == "ply") {
      Rvcg::vcgPlyWrite(mesh, filename = export_path, writeCol = FALSE)
    } else if (format == "obj") {
      Rvcg::vcgObjWrite(mesh, filename = export_path)
    } else if (format == "stl") {
      Rvcg::vcgSTLWrite(mesh, filename = export_path)
    }
    message("Mesh exported to: ", export_path)
  }
  
  return(mesh)
}


### Geometric measurement functions

## Compute depth metrics from aligned mesh
compute_depth <- function(mesh) {
  z_vals <- t(mesh$vb[3, ])
  abs(min(z_vals))
}


## Compute maximum distances along two principal axes in XY plane
compute_max_distance_axes <- function(mesh) {
  verts <- t(mesh$vb[1:3, ])
  xy_coords <- verts[, 1:2]  # Project to XY plane
  
  # Find maximum pairwise Euclidean distance
  dist_matrix <- as.matrix(dist(xy_coords))
  max_dist <- max(dist_matrix)
  
  # Find which two points define this maximum
  max_idx <- which(dist_matrix == max_dist, arr.ind = TRUE)[1,]
  p1 <- xy_coords[max_idx[1], ]
  p2 <- xy_coords[max_idx[2], ]
  
  # Vector defining longest axis
  longest_vec <- p2 - p1
  longest_direction <- longest_vec / sqrt(sum(longest_vec^2))
  
  # Perpendicular direction (rotate 90° in XY plane)
  perp_direction <- c(-longest_direction[2], longest_direction[1])
  
  # Project all points onto these directions and measure range
  proj_longest <- xy_coords %*% longest_direction
  proj_perp <- xy_coords %*% perp_direction
  
  axis_1 <- max(proj_longest) - min(proj_longest)
  axis_2 <- max(proj_perp) - min(proj_perp)
  
  return(data.frame(
    Max_Axis_1 = axis_1,
    Max_Axis_2 = axis_2
  ))
}

### Statistical analysis functions

## Robust Coefficient of Variation (rCV)
rCV <- function(x, na.rm = TRUE) {
  med <- median(x, na.rm = na.rm)
  mad_val <- mad(x, constant = 1, na.rm = na.rm)  # unscaled MAD
  rcv_val <- 100 * (1.4826 * mad_val / med)
  return(rcv_val)
}


# Import data -------------------------------------------------------------

## Configuration
wd <-  getwd()


# Volume calculations (pymeshlab)
vol <- read.csv(file.path(wd, "Volume.csv"))
vol_hull <- read.csv(file.path(wd, "Volume_hull.csv"))


# 3D meshes
surfaces_folder <- file.path(wd, "Components")
cropped_folder <- file.path(wd, "Cropped")
hulls_folder <- file.path(wd, "Capped_surfaces")


## Get file lists
surfaces_files <- list.files(path = surfaces_folder, pattern = "\\.ply$", full.names = TRUE)
cropped_files <- list.files(path = cropped_folder, pattern = "\\.ply$", full.names = TRUE)
hulls_files <- list.files(path = hulls_folder, pattern = "\\.ply$", full.names = TRUE)

if (length(surfaces_files) == 0) {
  stop("No PLY files found in Components folder")
}
if (length(cropped_files) == 0) {
  stop("No PLY files found in Cropped folder")
}
if (length(hulls_files) == 0) {
  stop("No PLY files found in Capped_surfaces folder")
}


## Load meshes
surfaces <- lapply(surfaces_files, function(file) {
  vcgPlyRead(file, updateNormals = TRUE)
})

cropped <- lapply(cropped_files, function(file) {
  vcgPlyRead(file, updateNormals = TRUE)
})

hulls <- lapply(hulls_files, function(file) {
  vcgPlyRead(file)
})



## Set names
names(surfaces) <- tools::file_path_sans_ext(basename(surfaces_files))
names_surfaces <- names(surfaces)

names(cropped) <- tools::file_path_sans_ext(basename(cropped_files))
names_cropped <- names(cropped)

names(hulls) <- tools::file_path_sans_ext(basename(hulls_files))
names_hulls <- names(hulls)


# Prepare dataframe -------------------------------------------------------

## Create data frame
metrics <- data.frame(Mesh_file = names_surfaces)


# Clean-up names
metrics$Name <- sub("_remeshed.*$", "", metrics$Mesh_file)
metrics$Name <- substring(metrics$Name, 4)


# Add site and specimen columns
metrics$Site   <- sub("^(...).*", "\\1", metrics$Name)
metrics$Specimen <- sub("^[A-Za-z]+_(.*)$", "\\1", metrics$Name)


# Add sample name (dev / val)
metrics$Sample <- ifelse(
  grepl("^(Bur|Mau|Scl)", metrics$Site), "Validation", "Development")


# Map site codes to full names
site_code <- substr(metrics$Site, 1, 3)
site_map <- c(Bur = "Burnot", 
              Scl = "Sclaigneaux", 
              Mau = "Maurenne", 
              Isl = "Isle Adam", 
              Pet = "Petit Morin")

metrics$Site <- site_map[site_code]


# Add side
metrics$Side <- "Right"


# Quantification of extracted surfaces ------------------------------------

### Align surfaces to boundary plane
dir.create(file.path(wd, "aligned_models"), showWarnings = FALSE, recursive = TRUE)

aligned <- mapply(
  function(mesh, name) {
    export_file <- file.path(wd, paste0("aligned_models/", name, "_aligned.ply"))
    align_to_boundary_plane(mesh, export_path = export_file)
  },
  mesh = surfaces,
  name = names(surfaces),
  SIMPLIFY = FALSE
)


### Size metrics

## Compute maximum depth
metrics$Depth <- sapply(aligned, compute_depth)
metrics$Depth <- as.numeric(metrics$Depth)


## Compute 2D footprint area
metrics$area2d <- sapply(aligned, area2d)


## Add volume previously computed with PyMeshLab
vol$Name <- substring(
  sub("_remeshed.*$", "", substr(vol$Mesh.File, 1, nchar(vol$Mesh.File) - 4)), 4)
metrics <- merge(metrics, vol[,c("Name","Volume")], by="Name", all.x = TRUE)


### Shape metrics

## Sphericity

# Calculate surface areas of hull
hull_areas <- sapply(hulls, vcgArea)


# Add hull areas to metrics
hull_metrics <- data.frame(Name = names(hull_areas), Hull_area = as.numeric(hull_areas))
hull_metrics$Name <- substring(
  sub("_remeshed.*$", "", substr(hull_metrics$Name, 1, nchar(hull_metrics$Name) - 4)), 8)
metrics <- merge(metrics, hull_metrics, by="Name", all.x = TRUE)


# Sphericity calculation
metrics$Sphericity <- (pi^(1/3) * (6 * metrics$Volume)^(2/3)) / metrics$Hull_area


## Axis ratios
max_dist_df <- lapply(aligned, compute_max_distance_axes)
metrics$first_axis_dist <- sapply(max_dist_df, function(x) x$Max_Axis_1)
metrics$second_axis_dist <- sapply(max_dist_df, function(x) x$Max_Axis_2)
metrics$Elongation <- metrics$second_axis_dist / metrics$first_axis_dist
metrics$Flatness <- metrics$Depth / metrics$first_axis_dist


### Surface metrics

## Solidity

# Hull volume
vol_hull$Name <- substring(
  sub("_remeshed.*$", "", substr(vol_hull$Mesh.File, 1, nchar(vol_hull$Mesh.File) - 4)), 4)
names(vol_hull)[names(vol_hull) == "Volume"] <- "Volume_hull"
metrics <- merge(metrics, vol_hull[,c("Name","Volume_hull")], by="Name", all.x = TRUE)


# Solidity calculation
metrics$Solidity <- metrics$Volume / metrics$Volume_hull


## Surface roughness (as median absolute deviation of mean curvature)

# Verify name matching
normalize_surface_name <- function(x) {
  x <- substring(x, 4)
  x <- sub("cropped.*$", "cropped", x)
  x
}

cbind(
  original = names(surfaces)[1:10],
  normalized = normalize_surface_name(names(surfaces)[1:10])
)

setdiff(normalize_surface_name(names(surfaces)), names(cropped))


# Compute mean curvature for full cropped meshes
full_curv <- lapply(cropped, function(mesh) {
  vcgCurve(mesh)$meanvb
})


# Extract full XYZ coordinates for nearest neighbor mapping
full_xyz <- lapply(cropped, function(mesh) {
  t(mesh$vb[1:3, ])
})


# Map curvature from cropped mesh to extracted surfaces
map_curvature <- function(lesion_mesh, full_xyz, full_curv) {
  lesion_xyz <- t(lesion_mesh$vb[1:3, ])
  nn <- get.knnx(full_xyz, lesion_xyz, k = 1)
  full_curv[nn$nn.index[, 1]]
}

surface_curv <- lapply(names(surfaces), function(name) {
  full_name <- normalize_surface_name(name)
  
  map_curvature(
    lesion_mesh = surfaces[[name]],
    full_xyz    = full_xyz[[full_name]],
    full_curv   = full_curv[[full_name]]
  )
})

names(surface_curv) <- names(surfaces)


# Calculate median absolute deviation of curvature for each surface
metrics$meanvb_mad <- sapply(surface_curv, function(curv) {
  if (all(is.na(curv))) return(NA)
  mad(curv, na.rm = TRUE)
})


# Descriptive statistics ------------------------------------------------------

## Exploration of variables

pairs <- GGally::ggpairs(metrics %>% 
                  select(Depth, area2d, Volume, Sphericity, Elongation, Flatness,
                         Solidity, meanvb_mad),
                aes(color = metrics$Sample, alpha = 0.7))

print(pairs)


## Descriptive statistics

# Summary
metrics_summary <- metrics %>%
  summarise(
    Count = n(),
    across(where(is.numeric), list(
      median = ~median(.x, na.rm = TRUE),
      Q1 = ~quantile(.x, 0.25, na.rm = TRUE),
      Q3 = ~quantile(.x, 0.75, na.rm = TRUE),
      IQR = ~IQR(.x, na.rm = TRUE),
      rCV = ~rCV(.x, na.rm = TRUE),
      min = ~min(.x, na.rm = TRUE),
      max = ~max(.x, na.rm = TRUE)
    ), .names = "{.col}_{.fn}")
  )


# Transpose for easier viewing
metrics_summary <- as.data.frame(t(metrics_summary))
metrics_summary <- rownames_to_column(metrics_summary, var = "Metric")
colnames(metrics_summary) <- c("Metric", "Value")
metrics_summary <- metrics_summary %>%
  mutate(Value = as.numeric(trimws(Value)))


## Raw data frame

# Clean-up for export
metrics$Name <- NULL
metrics$Mesh_file <- NULL


# Re-arrange columns
metrics <- metrics %>% select(Sample, Site, Specimen, everything())


# Reorder rows
sample_order <- c("Development", "Validation")

metrics <- metrics %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%
  arrange(Sample)
metrics <- metrics %>%
  group_by(Sample, Site) %>%
  arrange(as.numeric(gsub("[^0-9]", "", Specimen)), .by_group = TRUE) %>% 
  ungroup()


# Round values
metrics_rounded <- metrics %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
metrics_summary_rounded <- metrics_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


# Correlation between variables -------------------------------------------

## Select numeric metrics for correlation analysis
numeric_metrics <- metrics %>%
  select(Depth, area2d, Volume, Sphericity, Elongation, Flatness,
         Solidity, meanvb_mad) %>%
  rename(
    `2D Area` = area2d,
    `MAD of \nMean Curvature` = meanvb_mad
  )


## Compute Spearman correlation matrix
cor_res <- rcorr(as.matrix(numeric_metrics), type = "spearman")
cor_matrix <- cor_res$r
p_matrix <- cor_res$P

print(round(cor_res$r, 2))

                                   
## Visualize correlation matrix

# Color palette
pal <- vangogh_palette("SunflowersMunich")


# Add significance stars
sig_labels <- ifelse(p_matrix < 0.001, "***",
                     ifelse(p_matrix < 0.01, "**",
                            ifelse(p_matrix < 0.05, "*", "")))
sig_labels[is.na(sig_labels)] <- ""


# Remove diagonal
diag(cor_matrix) <- NA
diag(p_matrix) <- NA
diag(sig_labels) <- ""


# Mask triangles
cor_upper <- cor_matrix
cor_upper[lower.tri(cor_upper)] <- NA

p_lower <- p_matrix
p_lower[upper.tri(p_lower)] <- NA 


# Melt matrices
cor_melt <- melt(cor_upper, na.rm = TRUE)
p_melt <- melt(p_lower, na.rm = TRUE)
sig_melt <- melt(sig_labels)


# Combine for plotting
plot_data <- merge(cor_melt, p_melt, by = c("Var1", "Var2"), all = TRUE)
plot_data <- merge(plot_data, sig_melt, by = c("Var1", "Var2"), all = TRUE)
names(plot_data) <- c("Var1", "Var2", "cor", "pval", "sig")


# Replace NA stars by empty string
plot_data$sig[is.na(plot_data$sig)] <- ""


# Create labels
plot_data$label <- ifelse(
  !is.na(plot_data$cor),
  paste0(
    formatC(plot_data$cor, format = "f", digits = 2),
    plot_data$sig            # stars ONLY here
  ),
  ifelse(
    !is.na(plot_data$pval),
    ifelse(
      plot_data$pval < 0.001,
      "<0.001",              # NO stars for p-values
      formatC(plot_data$pval, format = "f", digits = 3)
    ),
    ""
  )
)


# Reverse y-axis order
plot_data$Var2 <- factor(
  plot_data$Var2,
  levels = rev(levels(plot_data$Var2))
)


# Plot
matrix_sig <- ggplot(plot_data, aes(Var1, Var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(low = pal[1], mid = "white", high = pal[5], midpoint = 0, 
                       limits = c(-1, 1), na.value = "grey80") +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Correlation\ncoefficient") +
  theme(legend.title = element_text(margin = margin(b = 12)))

ggsave("corrmatrix_sign.png", plot = matrix_sig, path = wd, width = 10, height = 6)


# Export results ----------------------------------------------------------

## Export data frames

wb <- createWorkbook()
addWorksheet(wb, "Metrics_full")
writeData(wb, "Metrics_full", metrics)
addWorksheet(wb, "Metrics_rounded")
writeData(wb, "Metrics_rounded", metrics_rounded)
addWorksheet(wb, "Metrics_summary")
writeData(wb, "Metrics_summary", metrics_summary)
addWorksheet(wb, "Metrics_summary_rounded")
writeData(wb, "Metrics_summary_rounded", metrics_summary_rounded)
addWorksheet(wb, "Cor_matrix")
writeData(wb, "Cor_matrix", cor_matrix)
addWorksheet(wb, "p_matrix")
writeData(wb, "p_matrix", p_matrix)
saveWorkbook(wb, file = file.path(wd, "Results_metrics.xlsx"), overwrite = TRUE)


