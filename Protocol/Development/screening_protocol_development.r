############################################################################
# Title: Microtrauma at the Humeral Medial Epicondyle - Screening Protocol Development
# 
# Author: Elle B. K. Liagre
# Email: elle.liagre@u-bordeaux.fr
# ORCID: https://orcid.org/0000-0002-8993-3266
# Date: 2025-12-18
# 
# Description: 
#   This script develops a two-stage classification protocol for identifying
#   microtrauma lesions at the humeral medial epicondyle using 3D mesh analysis.
#   The protocol combines spatial filtering (convex hull with buffer zone) and
#   threshold-based classification (depth metric) to distinguish true lesions
#   from other surface features.
#
# Requirements:
#    - R (>= 4.3.2)
#    - Rvcg (>= 0.25)
#    - dplyr (>= 1.1.4) 
#    - tidyr (>= 1.3.1)
#    - ggplot2 (>= 3.5.2)
#    - openxlsx (>= 4.2.5.2)
#    - pracma (>= 2.4.4)
#    - rgl (>= 1.3.18)
#    - doolkit (>= 1.42.2)
#    - geometry (>= 0.5.2)
#    - FNN (>= 1.1.4.1)
#    - mesheR (>= 0.4.200213)
#    - e1071 (>= 1.7.16)
#    - effsize (>= 0.8.1)
#
# Input data:
#     - Extracted surface components as PLY files
#     - Cropped humeri as PLY files
#
# License: This script is licensed under the GNU General Public License v3.0.
# See https://www.gnu.org/licenses/gpl-3.0.html for more details.
#
###########################################################################


# Clean workspace

rm(list=ls())


# Dependencies ------------------------------------------------------------

library(Rvcg)
library(openxlsx)
library(dplyr)
library(tidyr)
library(pracma)
library(rgl)
library(ggplot2)
library(doolkit)
library(geometry)
library(FNN)
library(mesheR)
library(e1071)
library(effsize)


# Functions ---------------------------------------------------------------

### Utility functions

## Extract clean specimen name from file name
extract_clean_name <- function(x) {
  # Remove "AO_" prefix
  x_no_prefix <- sub("^AO_", "", x)
  
  # Get everything before "_remeshed"
  prefix <- vapply(strsplit(x_no_prefix, "_remeshed", fixed = TRUE),
                   `[`, character(1), 1)
  
  # Get the last "_compXX"
  suffix <- sub(".*(_comp\\d{1,2}(?:_\\d+)?)$", "\\1", x)
  
  # Combine
  paste0(prefix, suffix)
}


## Get bounding box of mesh
get_mesh_bbox <- function(mesh) {
  verts <- t(mesh$vb[1:3, ])
  list(
    min = apply(verts, 2, min),
    max = apply(verts, 2, max)
  )
}


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


## Normalize mesh coordinates to [0,1] range
normalize_mesh <- function(mesh, specimen_name, surfaces_df) {
  
  # Get bounding box for this specimen
  bb <- surfaces_df[surfaces_df$Name == specimen_name, 
                    c("xmin","xmax","ymin","ymax","zmin","zmax")]
  
  verts <- t(mesh$vb[1:3, ])
  
  # Normalize to [0,1]
  verts_norm <- cbind(
    (verts[,1] - bb$xmin) / (bb$xmax - bb$xmin),
    (verts[,2] - bb$ymin) / (bb$ymax - bb$ymin),
    (verts[,3] - bb$zmin) / (bb$zmax - bb$zmin)
  )
  
  # Update mesh
  mesh$vb[1:3, ] <- t(verts_norm)
  return(mesh)
}


### Geometric measurement functions

## Extract geometric centroid from mesh
extract_centroid <- function(mesh) {
  verts <- t(mesh$vb[1:3, ])
  faces <- t(mesh$it)
  
  # Ensure meshes has faces
  if (nrow(faces) == 0) {
    warning("Mesh has zero faces.")
    return(rep(NA, 3))
  }
  
  # Extract triangle vertex coordinates
  A <- verts[faces[, 1], , drop = FALSE]
  B <- verts[faces[, 2], , drop = FALSE]
  C <- verts[faces[, 3], , drop = FALSE]
  
  # Triangle centroids
  tri_centroids <- (A + B + C) / 3
  
  # Triangle areas using cross product
  crossprod_vec <- B - A
  crossprod_vec2 <- C - A
  cp <- crossprod_vec[,2]*crossprod_vec2[,3] - crossprod_vec[,3]*crossprod_vec2[,2]
  cp2 <- crossprod_vec[,3]*crossprod_vec2[,1] - crossprod_vec[,1]*crossprod_vec2[,3]
  cp3 <- crossprod_vec[,1]*crossprod_vec2[,2] - crossprod_vec[,2]*crossprod_vec2[,1]
  tri_areas <- 0.5 * sqrt(cp^2 + cp2^2 + cp3^2)
  
  # Remove zero-area triangles
  keep <- tri_areas > 0
  if (!any(keep)) {
    warning("Mesh has only zero-area faces.")
    return(colMeans(verts))  # fallback to mean of vertices
  }
  
  tri_centroids <- tri_centroids[keep, , drop = FALSE]
  tri_areas <- tri_areas[keep]
  
  # Weighted centroid
  centroid <- colSums(tri_centroids * tri_areas) / sum(tri_areas)
  as.numeric(centroid)
}


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

## Find optimal classification threshold
find_optimal_thresholds <- function(df, metric, direction = "greater") {
  values <- df[[metric]]
  thresholds <- quantile(values, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE)
  
  best_youden <- -Inf
  best_f1 <- -Inf
  best_youden_stats <- NULL
  best_f1_stats <- NULL
  
  for (thresh in unique(thresholds)) {
    # Apply classification rule
    if (direction == "greater") {
      predicted <- ifelse(df[[metric]] >= thresh, 1, 0)
    } else {
      predicted <- ifelse(df[[metric]] <= thresh, 1, 0)
    }
    
    # Confusion matrix
    tp <- sum(predicted == 1 & df$Lesion == 1, na.rm = TRUE)
    fp <- sum(predicted == 1 & df$Lesion == 0, na.rm = TRUE)
    tn <- sum(predicted == 0 & df$Lesion == 0, na.rm = TRUE)
    fn <- sum(predicted == 0 & df$Lesion == 1, na.rm = TRUE)
    
    # Performance metrics
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    youden <- sensitivity + specificity - 1
    
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- sensitivity
    f1 <- if ((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0
    
    # Track best Youden
    if (youden > best_youden) {
      best_youden <- youden
      best_youden_stats <- list(
        threshold = thresh,
        youden = youden,
        sensitivity = sensitivity,
        specificity = specificity,
        precision = precision,
        recall = recall,
        f1 = f1,
        tp = tp, fp = fp, tn = tn, fn = fn
      )
    }
    
    # Track best F1
    if (f1 > best_f1) {
      best_f1 <- f1
      best_f1_stats <- list(
        threshold = thresh,
        youden = youden,
        sensitivity = sensitivity,
        specificity = specificity,
        precision = precision,
        recall = recall,
        f1 = f1,
        tp = tp, fp = fp, tn = tn, fn = fn
      )
    }
  }
  
  return(list(
    best_youden = best_youden_stats,
    best_f1 = best_f1_stats
  ))
}


# Import data -------------------------------------------------------------

## Configuration
wd <-  getwd()

surfaces_folder <- file.path(wd, "Components")
cropped_folder <- file.path(wd, "Cropped")


## Get file lists
surfaces_files <- list.files(path = surfaces_folder, pattern = "\\.ply$", full.names = TRUE)
cropped_files <- list.files(path = cropped_folder, pattern = "\\.ply$", full.names = TRUE)


## Load meshes
surfaces <- lapply(surfaces_files, function(file) {
  vcgPlyRead(file, updateNormals = TRUE)
})

cropped <- lapply(cropped_files, function(file) {
  vcgPlyRead(file, updateNormals = TRUE)
})


## Set names
names(surfaces) <- tools::file_path_sans_ext(basename(surfaces_files))
names_surfaces <- names(surfaces)

names(cropped) <- tools::file_path_sans_ext(basename(cropped_files))
names_cropped <- names(cropped)


# Data encoding -----------------------------------------------------------

## Create metadata data frame
surfaces_df <- data.frame(Name_file = names_surfaces)
surfaces_df$Name <- extract_clean_name(surfaces_df$Name_file)
names(surfaces) <- surfaces_df$Name
surfaces_df$Specimen <- sub("^(?:AO_)?(.*?)(_comp.*)$", "\\1", surfaces_df$Name_file)


## Code lesion status
lesion_names <- c(
  "Bur_224_comp1_1", "Isl_634_60_comp2_1", "Mau_97_comp1_2",
  "Mau_111_comp1_2", "Pet_21186_comp2_1", "Pet_21203_comp2_2",
  "Pet_21208_comp3_1", "Pet_21210_comp3_1", "Pet_21211_comp1_2",
  "Pet_21212_comp1_1", "Pet_21220_comp10_1", "Pet_21229_comp2_1",
  "Scl_3_comp2_1", "Scl_27_comp1_2"
)

surfaces_df$Lesion <- ifelse(surfaces_df$Name %in% lesion_names, 1, 0)
surfaces_df$Lesion <- as.factor(surfaces_df$Lesion)


# Spatial characteristics ----------------------------------------------

## Extract centroids
centroids <- t(sapply(surfaces, extract_centroid))
surfaces_df$centroid_x <- centroids[, 1]
surfaces_df$centroid_y <- centroids[, 2]
surfaces_df$centroid_z <- centroids[, 3]


## Get bounding boxes from cropped meshes
bbox <- lapply(cropped, get_mesh_bbox)

bbox_df <- do.call(rbind, lapply(names(bbox), function(n) {
  b <- bbox[[n]]
  data.frame(
    Specimen = n,
    xmin = b$min[1], ymin = b$min[2], zmin = b$min[3],
    xmax = b$max[1], ymax = b$max[2], zmax = b$max[3]
  )
}))

surfaces_df <- merge(surfaces_df, bbox_df, by = "Specimen", all.x = TRUE)


## Normalize centroid coordinates to [0,1] range
surfaces_df$centroid_x_norm <- (surfaces_df$centroid_x - surfaces_df$xmin) /
  (surfaces_df$xmax - surfaces_df$xmin)

surfaces_df$centroid_y_norm <- (surfaces_df$centroid_y - surfaces_df$ymin) /
  (surfaces_df$ymax - surfaces_df$ymin)

surfaces_df$centroid_z_norm <- (surfaces_df$centroid_z - surfaces_df$zmin) /
  (surfaces_df$zmax - surfaces_df$zmin)


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
surfaces_df$Depth <- sapply(aligned, compute_depth)
surfaces_df$Depth <- as.numeric(surfaces_df$Depth)


## Compute 2D footprint area
surfaces_df$area2d <- sapply(aligned, function(x) area2d(x))


## Compute axis distances
max_dist_df <- lapply(aligned, compute_max_distance_axes)
surfaces_df$first_axis_dist <- sapply(max_dist_df, function(x) x$Max_Axis_1)
surfaces_df$second_axis_dist <- sapply(max_dist_df, function(x) x$Max_Axis_2)
surfaces_df$Elongation <- surfaces_df$second_axis_dist / surfaces_df$first_axis_dist
surfaces_df$Flatness <- surfaces_df$Depth / surfaces_df$first_axis_dist


### Shape metrics

## Surface roughness (as standard deviation of mean curvature)

# Verify name matching
normalize_full_name <- function(x) {
  sub("_remeshed.*$", "", x)
}

normalize_surface_name <- function(x) {
  sub("_comp.*$", "", x)
}


# Compute mean curvature for full cropped meshes
full_curv <- lapply(cropped, function(mesh) {
  vcgCurve(mesh)$meanvb
})


# Extract full XYZ coordinates for nearest neighbor mapping
full_xyz <- lapply(cropped, function(mesh) {
  t(mesh$vb[1:3, ])
})

names(full_xyz)  <- normalize_full_name(names(full_xyz))
names(full_curv) <- normalize_full_name(names(full_curv))


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


# Calculate standard deviation of curvature for each surface
surfaces_df$meanvb_sd <- sapply(surface_curv, function(curv) {
  if (all(is.na(curv))) return(NA)
  sd(curv, na.rm = TRUE)
})


# Create subsamples -------------------------------------------------------

## Remove NA from dataframe
surfaces_df <- na.omit(surfaces_df)


## Flatten any list columns
surfaces_df <- data.frame(lapply(surfaces_df, function(x) {
  if (is.list(x)) {
    # Flatten the list to a numeric vector
    unlist(x)
  } else {
    x
  }
}))


## Factorise "lesion"
surfaces_df$Lesion <- as.factor(surfaces_df$Lesion)


## Assign to development or validation sample based on collection
surfaces_df$Sample <- ifelse(
  grepl("^(Bur|Mau|Scl)", surfaces_df$Name),
  "val",
  ifelse(grepl("^(Isl|Pet)", surfaces_df$Name), "dev", NA)
)


# Split data frames
train_df <- surfaces_df[surfaces_df$Sample == "dev",]
valid_df <- surfaces_df[surfaces_df$Sample == "val",]


# Split mesh lists
dev <- surfaces[surfaces_df$Sample == "dev"]
val <- surfaces[surfaces_df$Sample == "val"]


# Spatial filtering -------------------------------------------------------

## Exploration of development sample
summary <- train_df %>%
  group_by(Lesion) %>%
  summarise(
    mean_x = mean(centroid_x_norm),
    mean_y = mean(centroid_y_norm),
    mean_z = mean(centroid_z_norm),
    sd_x = sd(centroid_x_norm),
    sd_y = sd(centroid_y_norm),
    sd_z = sd(centroid_z_norm)
  )

print(summary)


## Create convex hull

# Get names of lesion surfaces
lesion_names <- train_df$Name[train_df$Lesion == 1]


# Extract the corresponding meshes from 'dev'
lesion_meshes <- dev[names(dev) %in% lesion_names]


# Normalize all lesion meshes
lesion_meshes_norm <- mapply(
  normalize_mesh,
  mesh = lesion_meshes,
  specimen_name = names(lesion_meshes),
  MoreArgs = list(surfaces_df = train_df),
  SIMPLIFY = FALSE
)


# Combine all vertices from lesion surfaces
all_verts <- do.call(rbind, lapply(lesion_meshes_norm, 
                                   function(mesh) t(mesh$vb[1:3, ])))


# Create 3D convex hull
hull <- convhulln(all_verts, options = "FA")
hull_mesh <- to.mesh3d(hull)


## Calculate buffer

# Calculate buffer distance using k-nearest neighbors
pos_mat <- as.matrix(train_df[train_df$Lesion == 1, 
                              c("centroid_x_norm", "centroid_y_norm", 
                                "centroid_z_norm")])


nn <- get.knn(pos_mat, k = 2)
nn_dist <- nn$nn.dist[,2]
buffer <- median(nn_dist)

print(paste("Selected buffer size:", round(buffer, 4)))


# Create buffered hull using mesh offset
hull_mesh <- vcgClean(hull_mesh, sel = 1)
hull_mesh <- vcgClean(hull_mesh, sel = 7)

buffer_mesh <- meshOffset(hull_mesh, offset = buffer)
buffer_mesh <- vcgClean(buffer_mesh, sel = 1)


# Convert buffered mesh to convex hull for inside/outside testing
buffer_mesh_pts <- t(buffer_mesh$vb[1:3, ])
buffer_mesh_hull <- convhulln(buffer_mesh_pts, options = "FA")
buffer_mesh_hull_mesh <- to.mesh3d(buffer_mesh_hull)


# Visualization of hull and buffer with points
open3d()
shade3d(hull_mesh, color = "red", alpha = 0.5)
shade3d(buffer_mesh_hull_mesh, color = "blue", alpha = 0.3)

# Points
points_mat <- as.matrix(train_df[, c("centroid_x_norm",
                                     "centroid_y_norm",
                                     "centroid_z_norm")])

# Color points by in_buffer
point_colors <- ifelse(train_df$Lesion == 1, "green", "orange")

# Add points to the 3D scene
points3d(points_mat, color = point_colors, size = 5)


## Evaluation of spatial filter on development sample

# Check if points are inside hull and buffer
points_mat <- as.matrix(train_df[, c("centroid_x_norm",
                                     "centroid_y_norm",
                                     "centroid_z_norm")])

train_df$in_hull <- inhulln(hull, points_mat)
train_df$in_buffer <- inhulln(buffer_mesh_hull, points_mat)


# Compute signed distances to hull and buffersurfaces
result_hull <- vcgClost(x = points_mat, mesh = hull_mesh, sign = TRUE)
result_buffer <- vcgClost(x = points_mat, mesh = buffer_mesh, sign = TRUE)

train_df$distances_hull <- result_hull$quality
train_df$distances_buffer <- result_buffer$quality


# Evaluation of performance
conf_matrix_spatial_dev <- table(train_df$Lesion, train_df$in_buffer)
sens <- mean(train_df$in_buffer[train_df$Lesion == 1])
fp_rate <- mean(train_df$in_buffer[train_df$Lesion == 0])

message("Spatial filter - Development sample:")
message("  Sensitivity: ", round(sens, 3))
message("  False positive rate: ", round(fp_rate, 3))
print(conf_matrix_spatial_dev)


## Evaluation of spatial filter on validation sample

# Check if points are inside hull and buffer
points_mat_valid <- as.matrix(valid_df[, c("centroid_x_norm",
                                           "centroid_y_norm",
                                           "centroid_z_norm")])

valid_df$in_hull <- inhulln(hull, points_mat_valid)
valid_df$in_buffer <- inhulln(buffer_mesh_hull, points_mat_valid)


# Compute signed distances to hull and buffer surfaces
result_valid_hull <- vcgClost(x = points_mat_valid, 
                              mesh = hull_mesh, sign = TRUE)
result_valid_buffer <- vcgClost(x = points_mat_valid, 
                                mesh = buffer_mesh, sign = TRUE)

valid_df$distances_hull <- result_valid_hull$quality
valid_df$distances_buffer <- result_valid_buffer$quality


# Evaluation of performance
conf_matrix_spatial_val <- table(valid_df$Lesion, valid_df$in_buffer)
sens_val <- mean(valid_df$in_buffer[valid_df$Lesion == 1])
spec_val <- mean(!valid_df$in_buffer[valid_df$Lesion == 0])

message("Spatial filter - Validation sample:")
message("  Sensitivity: ", round(sens_val, 3))
message("  Specificity: ", round(spec_val, 3))
print(conf_matrix_spatial_val)


# Threshold-based classification ------------------------------------------

## Define metrics
metric_vars <- c("Depth", "area2d", "Elongation", "Flatness", "meanvb_sd")

## Summary statistics by class
for (metric in metric_vars) {
  cat("---", metric, "---\n")
  cat("Training Positives:\n")
  print(summary(train_df[[metric]][train_df$Lesion == 1]))
  cat("\nTraining Negatives:\n")
  print(summary(train_df[[metric]][train_df$Lesion == 0]))
  cat("\nValidation Positives:\n")
  print(summary(valid_df[[metric]][valid_df$Lesion == 1]))
  cat("\nValidation Negatives:\n")
  print(summary(valid_df[[metric]][valid_df$Lesion == 0]))
  cat("\n")
}


## Distributions by metric

# Calculate skewness for each metric
skews <- sapply(metric_vars, function(m) skewness(train_df[[m]], na.rm = TRUE))

print(skews)


# Visualise metric distributions
train_long <- train_df %>%
  dplyr::select(Name, Lesion, all_of(metric_vars)) %>%
  pivot_longer(cols = all_of(metric_vars), names_to = "Metric", 
               values_to = "Value")

set.seed(9)

metrics_dev <- ggplot(train_long, aes(x = Lesion, y = Value, fill = Lesion)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.3) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("0" = "gray", "1" = "red")) +
  labs(title = "Development Sample: Metric Distributions by Class") +
  theme_minimal()

print(metrics_dev)


## Calculate effect sizes

train_df$Lesion <- factor(train_df$Lesion, levels = c("1","0"))

# Cohen's d (parametric)
cohens_d <- lapply(metric_vars, function(m) {
  cohen.d(train_df[[m]] ~ train_df$Lesion)
})


# Hedges G (small sample correction)
hedges_g <- lapply(metric_vars, function(m) {
  cohen.d(train_df[[m]] ~ train_df$Lesion, hedges.correction = TRUE)
})


# Glass's delta
glass_delta <- lapply(metric_vars, function(m) {
  cohen.d(train_df[[m]] ~ train_df$Lesion, pooled = FALSE)
})


# Cliff's delta (non-parametric)
cliffs_delta <- lapply(metric_vars, function(m) {
  cliff.delta(train_df[[m]] ~ train_df$Lesion)
})


# Compile effect sizes into a data frame
effect_size_df <- data.frame(
  Metric = metric_vars,
  cohens_d = sapply(cohens_d, function(x) x$estimate),
  cohens_d_lowCI = sapply(cohens_d, function(x) x$conf.int[1]),
  cohens_d_highCI = sapply(cohens_d, function(x) x$conf.int[2]),
  cohens_d_magnitude = sapply(cohens_d, function(x) x$magnitude),
  hedges_g = sapply(hedges_g, function(x) x$estimate),
  hedges_g_lowCI = sapply(hedges_g, function(x) x$conf.int[1]),
  hedges_g_highCI = sapply(hedges_g, function(x) x$conf.int[2]),
  hedges_g_magnitude = sapply(hedges_g, function(x) x$magnitude),
  glass_delta = sapply(glass_delta, function(x) x$estimate),
  glass_delta_lowCI = sapply(glass_delta, function(x) x$conf.int[1]),
  glass_delta_highCI = sapply(glass_delta, function(x) x$conf.int[2]),
  glass_delta_magnitude = sapply(glass_delta, function(x) x$magnitude),
  cliffs_delta = sapply(cliffs_delta, function(x) x$estimate),
  cliffs_delta_lowCI = sapply(cliffs_delta, function(x) x$conf.int[1]),
  cliffs_delta_highCI = sapply(cliffs_delta, function(x) x$conf.int[2]),
  cliffs_delta_magnitude = sapply(cliffs_delta, function(x) x$magnitude)
)

print(effect_size_df)


## Calculate metric correlations (Spearman for robustness)
cor_matrix <- cor(train_df[, c("Depth", "area2d", "meanvb_sd", 
                               "Elongation", "Flatness")], 
    method="spearman")

print(cor_matrix)


## Threshold-determination

# Determine direction for each metric based on effect size
metric_directions <- ifelse(sapply(cliffs_delta, function(x) x$estimate) > 0, 
                            "greater", "less")

# Find optimal thresholds
optimal_thresholds <- list()

for (i in seq_along(metric_vars)) {
  metric <- metric_vars[i]
  direction <- metric_directions[i]
  
  cat("\n", metric, "(direction:", direction, ")\n")
  
  opt <- find_optimal_thresholds(train_df, metric, direction)
  optimal_thresholds[[metric]] <- opt
  
  # Print best Youden threshold
  youden <- opt$best_youden
  cat("  Best Youden threshold:", round(youden$threshold, 2), "\n")
  cat("    Sensitivity:", round(youden$sensitivity, 3), 
      ", Specificity:", round(youden$specificity, 3),
      ", Precision:", round(youden$precision, 3), 
      ", Recall:", round(youden$recall, 3), 
      ", F1:", round(youden$f1, 3), "\n")
  
  # Print best F1 threshold
  f1 <- opt$best_f1
  cat("  Best F1 threshold:", round(f1$threshold, 2), "\n")
  cat("    Sensitivity:", round(f1$sensitivity, 3), 
      ", Specificity:", round(f1$specificity, 3),
      ", Precision:", round(f1$precision, 3), 
      ", Recall:", round(f1$recall, 3), 
      ", F1:", round(f1$f1, 3), "\n")
}


## Apply depth threshold (best performer)
train_df$Filter_depth_y <- ifelse(
  train_df$Depth > optimal_thresholds$Depth$best_youden$threshold, 1, 0)

valid_df$Filter_depth_y <- ifelse(
  valid_df$Depth > optimal_thresholds$Depth$best_youden$threshold, 1, 0)

conf_matrix_depth_dev <- table(train_df$Lesion, train_df$Filter_depth_y)
conf_matrix_depth_val <- table(valid_df$Lesion, valid_df$Filter_depth_y)

message("\nDepth filter - Development sample:")
print(conf_matrix_depth_dev)
message("\nDepth filter - Validation sample:")
print(conf_matrix_depth_val)


# Combining filters -------------------------------------------------------

## Apply both filters
train_df$Filtered <- ifelse(train_df$in_buffer == TRUE & train_df$Filter_depth_y == 1, "1", "0")
train_df$Filtered_label <- with(train_df, 
  ifelse(in_buffer & Filter_depth_y == 1 & Lesion == "1", "tp",
  ifelse(in_buffer & Filter_depth_y == 1 & Lesion == "0", "fp", "tn")))

valid_df$Filtered <- ifelse(valid_df$in_buffer == TRUE & valid_df$Filter_depth_y == 1, "1", "0")
valid_df$Filtered_label <- with(valid_df,
  ifelse(in_buffer & Filter_depth_y == 1 & Lesion == "1", "tp",
  ifelse(in_buffer & Filter_depth_y == 1 & Lesion == "0", "fp", "tn")))

conf_matrix_combined_dev <- table(train_df$Lesion, train_df$Filtered)
conf_matrix_combined_val <- table(valid_df$Lesion, valid_df$Filtered)

message("\nCombined filter - Development sample:")
print(conf_matrix_combined_dev)
message("\nCombined filter - Validation sample:")
print(conf_matrix_combined_val)

# Combine datasets for export
final_df <- rbind(train_df, valid_df)


# Export of results -------------------------------------------------------

## Export to Excel
wb <- createWorkbook()
addWorksheet(wb, "Spatial_dev")
writeData(wb, "Spatial_dev", conf_matrix_spatial_dev)
addWorksheet(wb, "Spatial_val")
writeData(wb, "Spatial_val", conf_matrix_spatial_val)
addWorksheet(wb, "Skewness")
writeData(wb, "Skewness", skews)
addWorksheet(wb, "Effect_sizes")
writeData(wb, "Effect_sizes", effect_size_df)
addWorksheet(wb, "Correlation")
writeData(wb, "Correlation", cor_matrix)
addWorksheet(wb, "Depth_dev")
writeData(wb, "Depth_dev", conf_matrix_depth_dev)
addWorksheet(wb, "Depth_val")
writeData(wb, "Depth_val", conf_matrix_depth_val)
addWorksheet(wb, "Combined_dev")
writeData(wb, "Combined_dev", conf_matrix_combined_dev)
addWorksheet(wb, "Combined_val")
writeData(wb, "Combined_val", conf_matrix_combined_val)
addWorksheet(wb, "raw_data")
writeData(wb, "raw_data", final_df)

saveWorkbook(wb, file = file.path(wd, "Results_classification.xlsx"), 
             overwrite = TRUE)


## Export buffer as an RDS file
saveRDS(buffer_mesh, file = "buffer_mesh.rds")
saveRDS(buffer_mesh_hull, file = "buffer_mesh_hull.rds")


# Visualisation: setup ----------------------------------------------------

## Common plot elements
plot_theme <- theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "right"
  )

color_scale <- scale_color_manual(
  name = "",
  values = c("fp" = "grey60", "tp" = "#D74F1D", "tn" = "grey60"),
  labels = c("fp" = "False positive", "tp" = "True positive",
            "tn" = "True negative")
)

shape_scale <- scale_shape_manual(
  name = "",
  values = c("fp" = 16, "tp" = 17, "tn" = 1),
  labels = c("fp" = "False positive", "tp" = "True positive",
            "tn" = "True negative")
)

fill_scale <- scale_fill_manual(
  name = "",
  values = c("fp" = "grey60", "tp" = "#D74F1D", "tn" = "white"),
  labels = c("fp" = "False positive", "tp" = "True positive",
            "tn" = "True negative")
)


# Set view matrix
view_matrix <- matrix(
  c(
     0.9312624, -0.1320341,  0.3395843, 0,
     0.3642671,  0.3572133, -0.8600628, 0,
    -0.007746428, 0.924643517, 0.380754977, 0,
     0, 0, 0, 1
  ),
  nrow = 4,
  byrow = FALSE
)


# Visualisation: spatial filtering (2D views) (exploratory) ---------------

## Extract 3D vertices from the meshes

# Original hull
hull_pts <- t(hull_mesh$vb[1:3, ])  # N x 3 matrix

# Buffered hull
buffer_pts <- t(buffer_mesh$vb[1:3, ])  # M x 3 matrix


## Project to 2D for plotting
proj_xy <- hull_pts[, 1:2]
proj_buf_xy <- buffer_pts[, 1:2]

proj_xz <- hull_pts[, c(1,3)]
proj_buf_xz <- buffer_pts[, c(1,3)]

proj_yz <- hull_pts[, 2:3]
proj_buf_yz <- buffer_pts[, 2:3]


## Compute 2D convex hulls for plotting

# XY plane
hull_xy_idx <- chull(proj_xy); hull_xy_idx <- c(hull_xy_idx, hull_xy_idx[1])
buffer_xy_idx <- chull(proj_buf_xy); buffer_xy_idx <- c(buffer_xy_idx, buffer_xy_idx[1])

hull_edges_xy <- data.frame(x = proj_xy[hull_xy_idx,1],
                            y = proj_xy[hull_xy_idx,2])
buffer_edges_xy <- data.frame(x = proj_buf_xy[buffer_xy_idx,1], 
                              y = proj_buf_xy[buffer_xy_idx,2])

# XZ plane
hull_xz_idx <- chull(proj_xz); hull_xz_idx <- c(hull_xz_idx, hull_xz_idx[1])
buffer_xz_idx <- chull(proj_buf_xz); buffer_xz_idx <- c(buffer_xz_idx, buffer_xz_idx[1])

hull_edges_xz <- data.frame(x = proj_xz[hull_xz_idx,1], 
                            y = proj_xz[hull_xz_idx,2])
buffer_edges_xz <- data.frame(x = proj_buf_xz[buffer_xz_idx,1], 
                              y = proj_buf_xz[buffer_xz_idx,2])

# YZ plane
hull_yz_idx <- chull(proj_yz); hull_yz_idx <- c(hull_yz_idx, hull_yz_idx[1])
buffer_yz_idx <- chull(proj_buf_yz); buffer_yz_idx <- c(buffer_yz_idx, buffer_yz_idx[1])

hull_edges_yz <- data.frame(x = proj_yz[hull_yz_idx,1], 
                            y = proj_yz[hull_yz_idx,2])
buffer_edges_yz <- data.frame(x = proj_buf_yz[buffer_yz_idx,1], 
                              y = proj_buf_yz[buffer_yz_idx,2])


## Development sample - 2D plots with spatial filter

# Create separate variables for buffer status and lesion status
train_df$Spatial_filtered <- ifelse(train_df$in_buffer == TRUE & 
                                      train_df$Lesion ==1, "tp", 
                                    ifelse(train_df$in_buffer == TRUE & 
                                             train_df$Lesion == 0, "fp", "tn"))

# Create combined factor for proper ordering
train_df$Spatial_filtered <- factor(train_df$Spatial_filtered, 
                                 levels = c("tp", "fp", "tn"))


# XY plot
xy_dev <- ggplot(train_df, aes(x = centroid_x_norm, y = centroid_y_norm)) +
  geom_point(aes(color = Spatial_filtered, 
                 shape = Spatial_filtered,
                 fill = Spatial_filtered),
             alpha = 0.8, size = 1.8) +
  geom_polygon(data = buffer_edges_xy, aes(x=x, y=y), 
               inherit.aes = FALSE,
               fill = "#90B6A4", alpha = 0.2, 
               color = "#90B6A4", linetype = "dashed", linewidth = 1) +
  color_scale +
  shape_scale +
  fill_scale +
  labs(x = "Normalized x-coordinate",
       y = "Normalized y-coordinate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  plot_theme

print(xy_dev)


# YZ plot
yz_dev <- ggplot(train_df, aes(x = centroid_y_norm, y = centroid_z_norm)) +
  geom_point(aes(color = Spatial_filtered, 
                 shape = Spatial_filtered,
                 fill = Spatial_filtered),
             alpha = 0.8, size = 1.8) +
  geom_polygon(data = buffer_edges_yz, aes(x=x, y=y), 
               inherit.aes = FALSE,
               fill = "#90B6A4", alpha = 0.2, 
               color = "#90B6A4", linetype = "dashed", linewidth = 1) +
  color_scale +
  shape_scale +
  fill_scale +
  labs(x = "Normalized y-coordinate",
       y = "Normalized z-coordinate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  plot_theme

print(yz_dev)


# XZ
xz_dev <- ggplot(train_df, aes(x = centroid_x_norm, y = centroid_z_norm)) +
  geom_point(aes(color = Spatial_filtered, 
                 shape = Spatial_filtered,
                 fill = Spatial_filtered),
             alpha = 0.8, size = 1.8) +
  geom_polygon(data = buffer_edges_xz, aes(x=x, y=y), 
               inherit.aes = FALSE,
               fill = "#90B6A4", alpha = 0.2, 
               color = "#90B6A4", linetype = "dashed", linewidth = 1) +
  color_scale +
  shape_scale +
  fill_scale +
  labs(x = "Normalized x-coordinate",
       y = "Normalized z-coordinate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  plot_theme

print(xz_dev)


## Validation sample - 2D plots with spatial filter

# Create separate variables for buffer status and lesion status
valid_df$Spatial_filtered <- ifelse(valid_df$in_buffer == TRUE & valid_df$Lesion ==1, "tp", 
                                    ifelse(valid_df$in_buffer == TRUE & 
                                             valid_df$Lesion == 0, "fp", "tn"))

# Create combined factor for proper ordering
valid_df$Spatial_filtered <- factor(valid_df$Spatial_filtered, 
                                    levels = c("tp", "fp", "tn"))

# XY plot
xy_val <- ggplot(valid_df, aes(x = centroid_x_norm, y = centroid_y_norm)) +
  geom_point(aes(color = Spatial_filtered, 
                 shape = Spatial_filtered,
                 fill = Spatial_filtered),
             alpha = 0.8, size = 1.8) +
  geom_polygon(data = buffer_edges_xy, aes(x=x, y=y), 
               inherit.aes = FALSE,
               fill = "#90B6A4", alpha = 0.2, 
               color = "#90B6A4", linetype = "dashed", linewidth = 1) +
  color_scale +
  shape_scale +
  fill_scale +
  labs(x = "Normalized x-coordinate",
       y = "Normalized y-coordinate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  plot_theme

print(xy_val)


# YZ plot
yz_val <- ggplot(valid_df, aes(x = centroid_y_norm, y = centroid_z_norm)) +
  geom_point(aes(color = Spatial_filtered, 
                 shape = Spatial_filtered,
                 fill = Spatial_filtered),
             alpha = 0.8, size = 1.8) +
  geom_polygon(data = buffer_edges_yz, aes(x=x, y=y), 
               inherit.aes = FALSE,
               fill = "#90B6A4", alpha = 0.2, 
               color = "#90B6A4", linetype = "dashed", linewidth = 1) +
  color_scale +
  shape_scale +
  fill_scale +
  labs(x = "Normalized y-coordinate",
       y = "Normalized z-coordinate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  plot_theme

print(yz_val)


# XZ
xz_val <- ggplot(valid_df, aes(x = centroid_x_norm, y = centroid_z_norm)) +
  geom_point(aes(color = Spatial_filtered, 
                 shape = Spatial_filtered,
                 fill = Spatial_filtered),
             alpha = 0.8, size = 1.8) +
  geom_polygon(data = buffer_edges_xz, aes(x=x, y=y), 
               inherit.aes = FALSE,
               fill = "#90B6A4", alpha = 0.2, 
               color = "#90B6A4", linetype = "dashed", linewidth = 1) +
  color_scale +
  shape_scale +
  fill_scale +
  labs(x = "Normalized x-coordinate",
       y = "Normalized z-coordinate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  plot_theme

print(xz_val)


# Visualisation: spatial filtering (3D views) -----------------------------

## Development sample

# Prepare the data
train_colors <- ifelse(train_df$Spatial_filtered == "tp", "#D74F1D",
                       ifelse(train_df$Spatial_filtered == "fp", "grey40", "white"))


# Open 3D window
open3d(windowRect = c(0, 0, 1200, 1200))
bg3d("white")
par3d(mouseMode = "trackball")


# Draw buffer hull
shade3d(buffer_mesh_hull_mesh, color = "#90B6A4", alpha = 0.3, 
        specular = "black", shininess = 10)


# Draw points
train_points <- as.matrix(train_df[, c("centroid_x_norm",
                                       "centroid_y_norm",
                                       "centroid_z_norm")])


# Non-lesions as spheres
non_lesion_idx <- train_df$Lesion == 0
non_lesion_points <- train_points[non_lesion_idx, ]
non_lesion_colors <- train_colors[non_lesion_idx]
spheres3d(non_lesion_points, radius = 0.01, color = non_lesion_colors)


# Lesions as octahedrons
lesion_idx <- train_df$Lesion == 1
lesion_points <- train_points[lesion_idx, ]
lesion_colors <- train_colors[lesion_idx]
for(i in 1:nrow(lesion_points)) {
  tetra <- tetrahedron3d()
  tetra <- scale3d(tetra, 0.01, 0.01, 0.01)
  tetra <- translate3d(tetra, lesion_points[i,1], 
                       lesion_points[i,2], 
                       lesion_points[i,3])
  shade3d(tetra, col = lesion_colors[i])
}


# Invisible points to set full range
spheres3d(c(0.98, 0.98, 0.98), radius = 0.01, color = "white", alpha = 0)


# Add axes and labels
axes3d(edges = "bbox", labels = TRUE, tick = TRUE)
title3d(xlab = "Normalized x-coordinate", 
        ylab = "Normalized y-coordinate", 
        zlab = "Normalized z-coordinate")


# Add legend
legend3d("topright", 
         legend = c("True positive", "False positive", "True negative"),
         pch = c(17, 16, 16),
         col = c("#D74F1D", "grey40", "grey80"),
         cex = 1.2,
         bg = "white")

par3d(userMatrix = view_matrix)
rgl.snapshot(file.path(wd, "dev_3D.png"), fmt = "png")


## Validation sample

# Prepare the data
valid_colors <- ifelse(valid_df$Spatial_filtered == "tp", "#D74F1D",
                       ifelse(valid_df$Spatial_filtered == "fp", "grey40", "white"))


# Open 3D window
open3d(windowRect = c(0, 0, 1200, 1200))
bg3d("white")
par3d(mouseMode = "trackball")


# Draw buffer hull
shade3d(buffer_mesh_hull_mesh, color = "#90B6A4", alpha = 0.3, 
        specular = "black", shininess = 10)


# Draw points
valid_points <- as.matrix(valid_df[, c("centroid_x_norm",
                                       "centroid_y_norm",
                                       "centroid_z_norm")])


# Non-lesions as spheres
non_lesion_idx <- valid_df$Lesion == 0
non_lesion_points <- valid_points[non_lesion_idx, ]
non_lesion_colors <- valid_colors[non_lesion_idx]
spheres3d(non_lesion_points, radius = 0.01, color = non_lesion_colors)


# Lesions as octahedrons
lesion_idx <- valid_df$Lesion == 1
lesion_points <- valid_points[lesion_idx, ]
lesion_colors <- valid_colors[lesion_idx]
for(i in 1:nrow(lesion_points)) {
  tetra <- tetrahedron3d()
  tetra <- scale3d(tetra, 0.01, 0.01, 0.01)
  tetra <- translate3d(tetra, lesion_points[i,1], 
                       lesion_points[i,2], 
                       lesion_points[i,3])
  shade3d(tetra, col = lesion_colors[i])
}


# Invisible points to set full range
spheres3d(c(0.98, 0.98, 0.98), radius = 0.01, color = "white", alpha = 0)


# Add axes and labels
axes3d(edges = "bbox", labels = TRUE, tick = TRUE)
title3d(xlab = "Normalized x-coordinate", 
        ylab = "Normalized y-coordinate", 
        zlab = "Normalized z-coordinate")


# Add legend
legend3d("topright", 
         legend = c("True positive", "False positive", "True negative"),
         pch = c(17, 16, 16),
         col = c("#D74F1D", "grey40", "grey80"),
         cex = 1.2,
         bg = "white")


par3d(userMatrix = view_matrix)
rgl.snapshot(file.path(wd, "val_3D.png"), fmt = "png")


## 3D views with anatomical context
# Example with Pet_21228

# Normalise cropped mesh
mesh <- normalize_mesh(cropped$Pet_21228_remeshed_oriented_cropped,
                       "Pet_21228_comp1_1", surfaces_df)


# Open 3D window
open3d(windowRect = c(0, 0, 1200, 1200))
bg3d("white")
shade3d(mesh, color = "grey", alpha = 0.95, 
        specular = "black", shininess = 10)
shade3d(buffer_mesh_hull_mesh, color = "#90B6A4", alpha = 0.75, 
        specular = "black", shininess = 10)


# XY view
view3d(theta = 0, phi = 0, zoom = 0.8)
rgl.snapshot(file.path(wd, "buffer_3D_xy.png"), fmt = "png")


# YZ view
view3d(theta = 90, phi = 0, zoom = 0.8)
rgl.snapshot(file.path(wd, "buffer_3D_yz.png"), fmt = "png")


# XZ view
view3d(theta = 0, phi = -75, zoom = 0.8)
rgl.snapshot(file.path(wd, "buffer_3D_xz.png"), fmt = "png")


# Extra angled view
par3d(userMatrix = view_matrix)
rgl.snapshot(file.path(wd, "buffer_3D_extra.png"), fmt = "png") 


# Visualisation: threshold classification ---------------------------------

## Development sample - basic depth plot
set.seed(9)

depth_dev <- ggplot(train_df, aes(x = 0, y = Depth, 
                                  color = factor(Lesion), 
                                  shape = factor(Lesion))) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_youden$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_f1$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("rect", ymin = optimal_thresholds$Depth$best_youden$threshold, 
           ymax = Inf,
           xmin = -Inf, xmax = Inf,
           alpha = 0.2, fill = "#90B6A4") +
  annotate("text", x = 5.15, y = optimal_thresholds$Depth$best_youden$threshold + 0.1,
           label = "Youden Index", hjust = 0, color = "#90B6A4", size = 3.5) +
  annotate("text", x = 6, y = optimal_thresholds$Depth$best_f1$threshold + 0.1,
           label = "F1-score", hjust = 0, color = "#90B6A4", size = 3.5) +
  geom_jitter(position = position_jitter(width = 5, height = 0), size = 2, alpha = 0.8) +
  labs(x = "", y = "Depth (mm)", color = "Lesion") +
  scale_color_manual(name = "", values = c("0" = "grey60", "1" = "#D74F1D"),
                     labels = c("0" = "No lesion", "1" = "Lesion")) +
  scale_shape_manual(name = "", values = c("0" = 16, "1" = 17),
                     labels = c("0" = "No lesion", "1" = "Lesion")) +
  scale_x_continuous(limits = c(-7, 7)) +
  scale_y_continuous(limits = c(-0.05, 3.75)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

print(depth_dev)


## Development sample - with classification labels
train_df$Depth_filtered <- with(train_df, 
                                ifelse(Lesion == 0 & Depth < optimal_thresholds$Depth$best_youden$threshold, 
                                              "tn",
                                              ifelse(Lesion == 0 & Depth >= optimal_thresholds$Depth$best_youden$threshold, 
                                                     "fp",
                                                     "tp")))
train_df$Depth_filtered <- factor(train_df$Depth_filtered, levels = c("tp","fp", "tn"))

set.seed(9)

depth_dev_class <- ggplot(train_df, aes(x = 0, y = Depth, 
                                        color = factor(Depth_filtered), 
                                        shape = factor(Depth_filtered))) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_youden$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_f1$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("rect", ymin = optimal_thresholds$Depth$best_youden$threshold, 
           ymax = Inf,
           xmin = -Inf, xmax = Inf,
           alpha = 0.2, fill = "#90B6A4") +
  annotate("text", x = 5.45, y = optimal_thresholds$Depth$best_youden$threshold + 0.1,
           label = "Youden Index", hjust = 0, color = "#90B6A4", size = 3.5) +
  annotate("text", x = 6.2, y = optimal_thresholds$Depth$best_f1$threshold + 0.1,
           label = "F1-score", hjust = 0, color = "#90B6A4", size = 3.5) +
  geom_jitter(position = position_jitter(width = 5, height = 0), size = 2, alpha = 0.8) +
  labs(x = "", y = "Depth (mm)", color = "Lesion") +
  scale_color_manual(name = "", values = c("tn" = "grey60", 
                                           "fp" = "grey60",
                                           "tp" = "#D74F1D"),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_shape_manual(name = "", values = c("tn" = 1, 
                                           "fp" = 16, 
                                           "tp" = 17),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_x_continuous(limits = c(-7, 7)) +
  scale_y_continuous(limits = c(-0.05, 3.75)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

print(depth_dev_class)

ggsave("depth_dev_class.png", plot = depth_dev_class, path = wd, 
       width = 8, height = 6, dpi = 600)


## Validation sample - basic depth plot
valid_df$Lesion <- factor(valid_df$Lesion, levels = c("1","0"))

set.seed(9)

depth_val <- ggplot(valid_df, aes(x = 0, y = Depth, 
                                  color = factor(Lesion), 
                                  shape = factor(Lesion))) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_youden$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("rect", ymin = optimal_thresholds$Depth$best_youden$threshold, 
           ymax = Inf,
           xmin = -Inf, xmax = Inf,
           alpha = 0.2, fill = "#90B6A4") +
  geom_hline(yintercept = optimal_thresholds$Depth$best_f1$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("text", x = 5.15, y = optimal_thresholds$Depth$best_youden$threshold + 0.1,
           label = "Youden Index", hjust = 0, color = "#90B6A4", size = 3.5) +
  annotate("text", x = 6, y = optimal_thresholds$Depth$best_f1$threshold + 0.1,
           label = "F1-score", hjust = 0, color = "#90B6A4", size = 3.5) +
  geom_jitter(position = position_jitter(width = 5, height = 0), size = 2, alpha = 0.8) +
  labs(x = "", y = "Depth (mm)", color = "Lesion") +
  scale_color_manual(name = "", values = c("0" = "grey60", "1" = "#D74F1D"),
                     labels = c("0" = "No lesion", "1" = "Lesion")) +
  scale_shape_manual(name = "", values = c("0" = 16, "1" = 17),
                     labels = c("0" = "No lesion", "1" = "Lesion")) +
  scale_x_continuous(limits = c(-7, 7)) +
  scale_y_continuous(limits = c(-0.05, 3.75)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

print(depth_val)


## Validation sample - with classification labels
valid_df$Depth_filtered <- with(valid_df, ifelse(Lesion == 0 & Depth < optimal_thresholds$Depth$best_youden$threshold, 
                                              "tn",
                                              ifelse(Lesion == 0 & Depth >= optimal_thresholds$Depth$best_youden$threshold, 
                                                     "fp",
                                                     "tp")))
valid_df$Depth_filtered <- factor(valid_df$Depth_filtered, levels = c("tp","fp", "tn"))

set.seed(9)

depth_val_class <- ggplot(valid_df, aes(x = 0, y = Depth, 
                                        color = factor(Depth_filtered), 
                                        shape = factor(Depth_filtered))) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_youden$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("rect", ymin = optimal_thresholds$Depth$best_youden$threshold, 
           ymax = Inf,
           xmin = -Inf, xmax = Inf,
           alpha = 0.2, fill = "#90B6A4") +
  geom_hline(yintercept = optimal_thresholds$Depth$best_f1$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("text", x = 5.45, y = optimal_thresholds$Depth$best_youden$threshold + 0.1,
           label = "Youden Index", hjust = 0, color = "#90B6A4", size = 3.5) +
  annotate("text", x = 6.2, y = optimal_thresholds$Depth$best_f1$threshold + 0.1,
           label = "F1-score", hjust = 0, color = "#90B6A4", size = 3.5) +
  geom_jitter(position = position_jitter(width = 5, height = 0), size = 2, alpha = 0.8) +
  labs(x = "", y = "Depth (mm)", color = "Lesion") +
  scale_color_manual(name = "", values = c("tn" = "grey60", 
                                           "fp" = "grey60",
                                           "tp" = "#D74F1D"),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_shape_manual(name = "", values = c("tn" = 1, 
                                           "fp" = 16, 
                                           "tp" = 17),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_x_continuous(limits = c(-7, 7)) +
  scale_y_continuous(limits = c(-0.05, 3.75)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

print(depth_val_class)

ggsave("depth_val_class.png", plot = depth_val_class, path = wd, 
       width = 8, height = 6, dpi = 600)


# Visualisation: combined classification ----------------------------------

## Development sample: combined spatial and depth filtering

train_df$Filtered_label <- factor(train_df$Filtered_label, 
                                  levels = c("tp", "fp", "tn"))
dev_class <- ggplot(train_df, aes(x = distances_buffer, y = Depth, 
                                  colour = factor(Filtered_label), 
                                  shape = factor(Filtered_label))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_youden$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = 0, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("rect", ymin = optimal_thresholds$Depth$best_youden$threshold, 
           ymax = Inf,
           xmin = 0, xmax = Inf,
           alpha = 0.2, fill = "#90B6A4") +
  labs(x = "Distance to ROI (normalised units)", y = "Depth (mm)", 
       color = "Lesion") +
  scale_color_manual(name = "", values = c("tn" = "grey60", 
                                           "fp" = "grey60",
                                           "tp" = "#D74F1D"),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_shape_manual(name = "", values = c("tn" = 1, 
                                           "fp" = 16, 
                                           "tp" = 17),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_x_continuous(limits = c(-0.6, 0.2)) +
  scale_y_continuous(limits = c(-0.05, 3.75)) +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

print(dev_class)

ggsave("dev_class.png", plot = dev_class, path = wd, 
       width = 8, height = 6, dpi = 600)


## Validation sample: combined spatial and depth filtering

valid_df$Filtered_label <- factor(valid_df$Filtered_label, 
                                  levels = c("tp", "fp", "tn"))
val_class <- ggplot(valid_df, aes(x = distances_buffer, y = Depth, 
                                  colour = factor(Filtered_label),
                                  shape = factor(Filtered_label))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = optimal_thresholds$Depth$best_youden$threshold, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = 0, 
             color = "#90B6A4", linetype = "dashed", size = 1.2) +
  annotate("rect", ymin = optimal_thresholds$Depth$best_youden$threshold, 
           ymax = Inf,
           xmin = 0, xmax = Inf,
           alpha = 0.2, fill = "#90B6A4") +
  labs(x = "Distance to ROI (normalised units)", y = "Depth (mm)", 
       color = "Lesion") +
  scale_color_manual(name = "", values = c("tn" = "grey60", 
                                           "fp" = "grey60",
                                           "tp" = "#D74F1D"),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_shape_manual(name = "", values = c("tn" = 1, 
                                           "fp" = 16, 
                                           "tp" = 17),
                     labels = c("tn" = "True negative", 
                                "fp" = "False positive",
                                "tp" = "True positive")) +
  scale_x_continuous(limits = c(-0.6, 0.2)) +
  scale_y_continuous(limits = c(-0.05, 3.75)) +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

ggsave("val_class.png", plot = val_class, path = wd, 
       width = 8, height = 6, dpi = 600)

print(val_class)

