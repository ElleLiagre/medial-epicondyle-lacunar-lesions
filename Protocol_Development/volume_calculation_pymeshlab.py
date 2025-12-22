"""
Title: Microtrauma at the Humeral Medial Epicondyle - Volume Calculation Script
 
Author: Elle B. K. Liagre
Email: elle.liagre@u-bordeaux.fr
ORCID: https://orcid.org/0000-0002-8993-3266
Date: 2025-12-18

Description: 
    Batch processes 3D mesh files (.ply) to calculate volumes. The script
    verifies meshes are watertight and cleans them up. Cleaned models and
    a CSV file with volume measurements are saved in an 'Output' subfolder.
    
    More information can be found on:
    https://github.com/ElleLiagre/medial-epicondyle-lacunar-lesions

Requirements: 
    - Python 3.13+
    - pymeshlab (>= 2023.12.post3) (https://doi.org/10.5281/zenodo.4438750)

    Install with: pip install pymeshlab

Input data:
    - 3D mesh files in PLY format

License:
    GNU General Public License v3.0
    https://www.gnu.org/licenses/gpl-3.0.html
"""

# Import necessary libraries
import pymeshlab
import csv
import os

# ============================================================================
# CONFIGURATION - Update this path to your mesh folder
# ============================================================================

input_folder = r"path\to\input\folder"

# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Set up output folder
output_folder = os.path.join(input_folder, "Output")
output_csv = os.path.join(output_folder, "Volume.csv")
os.makedirs(output_folder, exist_ok=True)

# Get all .ply files
ply_files = [f for f in os.listdir(input_folder) if f.endswith(".ply")]

# Process each file and write results to CSV
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Mesh File", "Volume"])

    # Loop through each .ply file in the folder
    for ply_file in ply_files:
        try:
            mesh_path = os.path.join(input_folder, ply_file)
            output_mesh_path = os.path.join(output_folder, f"Vol_{ply_file}")

            # Load mesh
            ms = pymeshlab.MeshSet()
            ms.load_new_mesh(mesh_path)

            # Make mesh watertight
            ms.meshing_repair_non_manifold_vertices()
            ms.compute_selection_by_non_manifold_edges_per_face()
            ms.apply_selection_dilatation()
            ms.meshing_remove_selected_vertices_and_faces()
            ms.set_selection_none()

            ms.meshing_repair_non_manifold_vertices()
            ms.compute_selection_by_non_manifold_edges_per_face()
            ms.meshing_remove_selected_vertices_and_faces()
            ms.set_selection_none()

            ms.meshing_close_holes(maxholesize=3000, newfaceselected=False, refinehole=False)

            # Calculate volume (with retry logic if needed)
            max_attempts = 10
            attempts = 0
            volume = None

            while volume is None and attempts < max_attempts:
                vol = ms.get_geometric_measures()
                if 'mesh_volume' in vol:
                    volume = vol['mesh_volume']
                else:
                    # Additional cleanup and retry
                    ms.compute_selection_from_mesh_border()
                    ms.apply_selection_dilatation()
                    ms.meshing_remove_selected_vertices_and_faces()
                    ms.set_selection_none()
                    ms.meshing_close_holes(maxholesize=3000, newfaceselected=False, refinehole=False)

                attempts += 1

            if volume is None:
                print(f"Failed to compute volume for {ply_file} after 10 attempts.")
            else:
                print(f"Volume for {ply_file}: {volume} (after {attempts} attempts)")

                # Save processed mesh and write to CSV
                ms.save_current_mesh(output_mesh_path)
                writer.writerow([ply_file, volume])
                print(f"Processed and saved: {ply_file}")

        except Exception as e:
            print(f"Error processing {ply_file}: {e}")

print("Batch processing complete.")

