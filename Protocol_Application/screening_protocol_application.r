############################################################################
# Title: Microtrauma at the Humeral Medial Epicondyle - Screening Protocol Application
#
# Author: Elle B. K. Liagre
# Email: elle.liagre@u-bordeaux.fr
# ORCID: https://orcid.org/0000-0002-8993-3266
# Date: 2025-12-18
#
# Description:
#    This script applies a validated two-stage screening protocol to identify
#    microtrauma lesions on humeral medial epicondyle 3D surface meshes.
#    
#    The protocol combines:
#      1. Spatial filtering: Tests whether surface centroids fall within a 
#         predefined anatomical region of interest (buffered convex hull)
#      2. Depth-based classification: Applies a 0.58mm depth threshold to
#         distinguish true lesions from other surface features
#    
#    Surfaces passing both filters are classified as potential microtrauma lesions.
#    
#    Protocol development and validation are documented in the companion script
#    'screening_protocol_dev.R'. The pre-computed spatial filter (buffer_mesh_hull)
#    is loaded directly from the GitHub repository.
#
#    This script is designed for application to new datasets after the protocol
#    has been established and validated.
#
# Requirements:
#    - R (>= 4.3.2)
#    - Rvcg (>= 0.25)
#    - openxlsx (>= 4.2.5.2)
#    - geometry (geometry (>= 0.5.2))
#    - pracma (>= 2.4.4)
#
# Input data:
#    - Extracted surfaces in PLY format to screen
#    - Cropped humeral medial epicondyle meshes in PLY format
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
library(geometry)
library(pracma)


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


# Import data -------------------------------------------------------------

## Configuration
wd <- getwd()

surfaces_folder <- file.path(wd, "Components")
cropped_folder <- file.path(wd, "Cropped")


## Get file lists
surfaces_files <- list.files(path = surfaces_folder, pattern = "\\.ply$", full.names = TRUE)
cropped_files <- list.files(path = cropped_folder, pattern = "\\.ply$", full.names = TRUE)


if (length(surfaces_files) == 0) {
  stop("No PLY files found in Components folder")
}
if (length(cropped_files) == 0) {
  stop("No PLY files found in Cropped folder")
}


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


## Load pre-computed buffer mesh from GitHub repository
buffer <- readRDS(url("https://raw.githubusercontent.com/ElleLiagre/medial-epicondyle-lacunar-lesions/1c16b171d54e1595403e6782ad663bf6ad45548d/Protocol_Application/buffer_mesh_hull.rds"))


## Depth threshold
depth_threshold <- 0.58


# Data encoding -----------------------------------------------------------

## Create metadata data frame
df <- data.frame(Name_file = names_surfaces)
df$Name <- extract_clean_name(df$Name_file)
names(surfaces) <- df$Name
df$Specimen <- sub("^(?:AO_)?(.*?)(_comp.*)$", "\\1", df$Name_file)


# Spatial characteristics -------------------------------------------------

## Extract centroids
centroids <- t(sapply(surfaces, extract_centroid))
df$centroid_x <- centroids[, 1]
df$centroid_y <- centroids[, 2]
df$centroid_z <- centroids[, 3]


## Get bounding box from cropped meshes
bbox <- lapply(cropped, get_mesh_bbox)

bbox_df <- do.call(rbind, lapply(names(bbox), function(n) {
  b <- bbox[[n]]
  data.frame(
    Specimen = n,
    xmin = b$min[1], ymin = b$min[2], zmin = b$min[3],
    xmax = b$max[1], ymax = b$max[2], zmax = b$max[3]
  )
}))

df <- merge(df, bbox_df, by = "Specimen", all.x = TRUE)


## Normalize centroid coordinates to [0,1] range
df$centroid_x_norm <- (df$centroid_x - df$xmin) /
  (df$xmax - df$xmin)

df$centroid_y_norm <- (df$centroid_y - df$ymin) /
  (df$ymax - df$ymin)

df$centroid_z_norm <- (df$centroid_z - df$zmin) /
  (df$zmax - df$zmin)


# Depth computation -------------------------------------------------------

## Align surfaces to boundary plane
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


## Compute depth metrics
df$Depth <- sapply(aligned, compute_depth)
df$Depth <- as.numeric(df$Depth)


# Screening ---------------------------------------------------------------

## Spatial filtering

# Extract centroid points matrix
points_mat <- as.matrix(df[, c("centroid_x_norm", 
                               "centroid_y_norm", 
                               "centroid_z_norm")])


# Check if points are inside buffer
df$in_buffer <- inhulln(buffer, points_mat)



## Depth filtering
df$Filter_depth <- ifelse(df$Depth > depth_threshold, 1, 0)


## Combined filtering
df$Screened <- ifelse(df$in_buffer == TRUE & df$Filter_depth == 1, 1, 0)


# Export results ----------------------------------------------------------

## Save results to Excel
output_file <- file.path(wd, "screening_results.xlsx")
write.xlsx(df, file = output_file, rowNames = FALSE)


