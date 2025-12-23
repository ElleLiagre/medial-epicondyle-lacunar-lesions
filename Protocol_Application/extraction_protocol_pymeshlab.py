"""
Title: Microtrauma at the Humeral Medial Epicondyle - Combined Extraction Protocol
 
Author: Elle B. K. Liagre
Email: elle.liagre@u-bordeaux.fr
ORCID: https://orcid.org/0000-0002-8993-3266
Date: 2025-12-22

Description: 
    Combined pipeline that batch processes 3D mesh files (.ply) through:
    1. Ambient occlusion-based surface extraction
    2. Curvature-based refinement on extracted components

    More information can be found on:
    https://github.com/ElleLiagre/medial-epicondyle-lacunar-lesions

Requirements: 
    - Python ( >= 3.13.2)
    - pymeshlab (>= 2023.12.post3) (https://doi.org/10.5281/zenodo.4438750)
    - psutil (>= 7.1.0)

    Install with: pip install pymeshlab psutil

Input data:
    - One or more triangular surface meshes (.ply format) 
      representing cropped medial epicondyles of the humerus.

License:
    GNU General Public License v3.0
    https://www.gnu.org/licenses/gpl-3.0.html
"""

# Import necessary libraries
import pymeshlab
import os
import gc
import psutil
import time


# ============================================================================
# CONFIGURATION - Update this path to your mesh folder
# ============================================================================

input_folder = r"path\to\input\folder"

if not os.path.isdir(input_folder):
    raise ValueError(
        f"Input folder does not exist: {input_folder}\n"
        "Please update 'input_folder' in the configuration section."
    )


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def log_memory():
    """Log current memory usage"""
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024 ** 2
    print(f"[Memory usage: {mem_mb:.1f} MB]")


def process_ambient_occlusion(mesh_path, output_folder, ply_file):
    """
    Ambient occlusion-based surface extraction
    Returns: list of output component file paths
    """
    print(f"\n=== Ambient Occlusion Processing - {ply_file} ===")
    
    output_components = []
    
    try:
        # Create a new MeshSet
        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(mesh_path)

        # Step 1: Compute ambient occlusion
        ms.compute_scalar_ambient_occlusion_gpu(occmode='per-Vertex')

        # Step 2: Smooth color using Laplacian smoothing
        ms.apply_color_laplacian_smoothing_per_vertex(iteration=15)
        
        # Step 3: Adjust levels
        ms.apply_color_level_adjustment_per_vertex(in_min=165, in_max=0)

        # Step 4A: Select faces by color and extract components
        ms.compute_selection_by_color_per_face() 
        ms.apply_selection_inverse(invfaces=True, invverts=True)
        ms.meshing_remove_selected_vertices_and_faces()

        if ms.current_mesh() is None or ms.current_mesh().vertex_number() == 0 or ms.current_mesh().face_number() == 0:
            print(f"Skipping {ply_file} — no region matched selection criteria.")
            return output_components

        # Smooth borders
        ms.compute_selection_from_mesh_border()
        ms.apply_coord_laplacian_smoothing(selected=True)

        # Extract connected components from selection
        ms.generate_splitting_by_connected_components(delete_source_mesh=True)
        component_count = ms.mesh_number()

        if component_count > 0:
            print(f"{component_count} component(s) found in selection for {ply_file}")

        for i in range(1, component_count + 1):
            try:
                ms.set_current_mesh(i)

                # Save the processed mesh
                comp_filename = f"AO_{os.path.splitext(ply_file)[0]}_comp{i}.ply"
                comp_path = os.path.join(output_folder, comp_filename)
                ms.save_current_mesh(comp_path)
                output_components.append(comp_path)
                
                print(f"Saved component {i}: {comp_filename}")

            except Exception as comp_err:
                print(f"Error processing component {i} of {ply_file}: {comp_err}")
                
    except Exception as e:
        print(f"Error in ambient occlusion processing {ply_file}: {e}")
    
    finally:
        if "ms" in locals():
            del ms
        gc.collect()

    
    return output_components


def process_curvature(mesh_path, output_folder, component_file, max_iterations=20):
    """
    Curvature-based refinement
    """
    print(f"\n  → Curvature Processing - {os.path.basename(component_file)}")
    
    try:
        # Create a new MeshSet
        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(mesh_path)

        # Clean up mesh
        ms.compute_normal_per_vertex()
        ms.compute_normal_per_face()
        ms.meshing_re_orient_faces_coherently()

        if ms.current_mesh().vertex_number() < 10 or ms.current_mesh().face_number() < 5:
            print(f"  Skipping — mesh too small for curvature analysis.")
            return

        # Step 4B: Curvature-based selection
        ms.compute_curvature_and_color_apss_per_vertex(
            filterscale=10,
            sphericalparameter=5,
            curvaturetype='K1'
        )
        ms.compute_selection_by_scalar_per_vertex(minq=0.85, maxq=2.0)

        # Step 5: Refinement of delineation
        # Mark selection as quality 4, others as 0
        ms.compute_scalar_by_function_per_vertex(q='vsel ? 4 : 0')

        # Select border
        ms.set_selection_none(allfaces=True, allverts=True)
        ms.compute_selection_from_mesh_border()

        # Add 1 to border vertex quality
        ms.compute_scalar_by_function_per_vertex(q='vsel ? q + 1 : q')

        # Transfer quality to faces
        ms.compute_scalar_transfer_vertex_to_face()
        ms.set_selection_none(allfaces=True, allverts=True)

        EPS = 1e-6  # small tolerance for float comparison

        # Flood-fill iterations
        for i in range(max_iterations):
            faces_before = ms.current_mesh().selected_face_number()
            ms.set_selection_none(allfaces=True, allverts=True)
            
            # Select faces with quality in target range
            ms.compute_selection_by_scalar_per_face(minq=4.01 - EPS, maxq=5.0 + EPS)
            
            selected_faces = ms.current_mesh().selected_face_number()
            print(f"  Iteration {i+1}: selected faces = {selected_faces}")

            if selected_faces == 0:
                print("  No faces selected, skipping remaining flood-fill iterations.")
                break

            # Dilate selection to expand region
            ms.apply_selection_dilatation()

            # Mark dilated faces with quality 5
            ms.compute_scalar_by_function_per_face(q='fsel && fq > 1 ? 5 : fq')
            
            faces_after = ms.current_mesh().selected_face_number()

            if faces_before == faces_after:
                print(f"  Flooding complete after {i+1} iterations")
                break
            
            if i == max_iterations - 1:
                print(f"  Reached max iterations")

        # Remove high-quality faces
        ms.compute_selection_by_scalar_per_face(minq=4.01 - EPS, maxq=5.0 + EPS)
        selected_faces = ms.current_mesh().selected_face_number()

        if selected_faces > 0:
            print(f"  Removing {selected_faces} high-quality faces")
            ms.meshing_remove_selected_vertices_and_faces()
        else:
            print("  No high-quality faces found — skipping removal")

        # Check if mesh is empty
        if (ms.current_mesh() is None or 
            ms.current_mesh().vertex_number() == 0 or 
            ms.current_mesh().face_number() == 0):
            print(f"  No region matched selection criteria after curvature filtering.")
            return

        ms.compute_selection_from_mesh_border()
        ms.apply_coord_laplacian_smoothing(selected=True)

        # Extract connected components
        ms.generate_splitting_by_connected_components(delete_source_mesh=True)
        component_count = ms.mesh_number()
        print(f"  {component_count} final component(s) found")

        # Process each component
        for i in range(1, component_count + 1):
            try:
                ms.set_current_mesh(i)
                comp_filename = f"{os.path.splitext(os.path.basename(component_file))[0]}_curv{i}.ply"
                comp_path = os.path.join(output_folder, comp_filename)
                ms.save_current_mesh(comp_path)
                print(f"  Saved: {comp_filename}")

            except Exception as comp_err:
                print(f"  Error saving component {i}: {comp_err}")

    except Exception as e:
        print(f"  Error in curvature processing: {e}")

    finally:
        if "ms" in locals():
            del ms
        gc.collect()



# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    print("="*80)
    print("COMBINED EXTRACTION PROTOCOL")
    print("="*80)
    
    # Set up output folders
    ao_folder = os.path.join(input_folder, "AO")
    curv_folder = os.path.join(input_folder, "Curv")
    os.makedirs(ao_folder, exist_ok=True)
    os.makedirs(curv_folder, exist_ok=True)

    # Get list of input .ply files
    ply_files = sorted(f for f in os.listdir(input_folder) if f.endswith(".ply"))
    print(f"\nFound {len(ply_files)} .ply files to process.\n")

    if len(ply_files) == 0:
        print("No .ply files found in input folder!")
        return

    # Process each file through the pipeline
    for index, ply_file in enumerate(ply_files, start=1):
        print(f"\n{'='*80}")
        print(f"[{index}/{len(ply_files)}] PROCESSING: {ply_file}")
        print(f"{'='*80}")
        
        mesh_path = os.path.join(input_folder, ply_file)
        
        # Ambient Occlusion Processing
        ao_components = process_ambient_occlusion(mesh_path, ao_folder, ply_file)
        
        if len(ao_components) == 0:
            print(f"\nNo components extracted from Part 1 for {ply_file}")
            continue
        
        print(f"\n=== Processing {len(ao_components)} AO component(s) ===")
        
        # Process each AO component through curvature analysis
        for ao_comp_path in ao_components:
            process_curvature(ao_comp_path, curv_folder, ao_comp_path)
            time.sleep(0.1)  # Brief pause between components
        
        # Memory cleanup between files
        gc.collect()
        log_memory()
        time.sleep(0.3)
    
    print("\n" + "="*80)
    print("BATCH PROCESSING COMPLETE")
    print("="*80)
    print(f"AO components: {ao_folder}")
    print(f"Curvature refined components: {curv_folder}")


if __name__ == "__main__":
    main()