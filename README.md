# Microtrauma at the Humeral Medial Epicondyle: Quantifying Lacunar Lesions

This repository contains all scripts and protocols used in the study "Microtrauma at the Humeral Medial Epicondyle: Quantifying Lacunar Lesions to Reconstruct Activity Patterns in Past Populations."

<br> 
<br>

**Author and contact details**  
Author &nbsp;&nbsp;&nbsp; *Elle B. K. Liagre*  
Email &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *elle.liagre@u-bordeaux.fr*  
ORCID <span>&nbsp;&nbsp;&nbsp;</span>   *[https://orcid.org/0000-0002-8993-3266](https://orcid.org/0000-0002-8993-3266)*

<br>

---

## Table of Contents

- [Description](#description)
- [Dependencies and libraries](#dependencies-and-libraries)
- [How to Use](#how-to-use)
- [Data Availability](#data-availability)
- [Contributing](#contributing)
- [License and Citation](#license-and-citation)
<br>

---

## Description

This project introduces a semi-automated, standardized 3D cropping protocol for analyzing the medial epicondyle of the human humerus, with a focus on its entheses. The protocol, implemented in **3DSlicer** using the **PyMeshLab** package and the **SlicerMorph** extensions, standardizes input models —such as mesh resolution and anatomical orientation—ensuring high repeatability and reproducibility in the analysis. Anatomical landmarks define the boundaries of the medial epicondyle, accounting for variations in humeral size, and ensuring consistent capture of the entire entheseal surface, even amidst anatomical variability.

To improve efficiency, batch processing is incorporated, enabling the protocol to be applied to multiple humeral specimens in a streamlined manner. This feature is particularly beneficial for handling larger datasets, ensuring consistent results across specimens while minimizing manual intervention. The methodology enhances the accuracy and comparability of entheseal surface analysis, offering a valuable tool for studies on skeletal adaptations, physical activity, and human evolution, with strong potential for large-scale applications.

This protocol is presented in the study by **Liagre EBK, Remy F, Villotte S, Knüsel CJ. A Standardized, Three-Dimensional Cropping Protocol for Analyzing the Medial Epicondyle of the Humerus. Am J Biol Anthropol. 2025 Aug;187(4):e70100. doi: 10.1002/ajpa.70100.** Methodological choices and further information can be found there.
<br>
<br>
<br>

## Repository Structure
```
├── Protocol/
│   ├── Application/
│   │   ├── extraction_application.py             # Python script for 3D surface extraction
│   │   ├── screening_protocol_application.R      # R script for lesion screening
│   │   └── buffer_mesh_hull.RDS                  # 
│   └── Development
│       └── screening_protocol_development.R      # R script for developing the extraction protocol
├── Analysis/
│   ├── lesion_analysis.R                         # Statistical analysis of lesion morphology
│   ├── reliability_study.R                       # Intra-observer reliability analysis
│   └── volume_calculation.py                     # Python script for volume measurements
└── README.md
```


## Dependencies and libraries
This project was built with the following software and library versions. Please ensure these dependencies are installed to maintain compatibility and ensure proper functionality of the protocol. Refer to the [How to Use](#how-to-use) section for installation steps.

### Software and Extensions

### Python Libraries
- **PyMeshLab**: >= 2023.12.post2

### Citations

- Muntoni, A., & Cignoni, P. (2024). *PyMeshLab: PyMeshLab v2023.12.post2*. Zenodo. [doi: 10.5281/zenodo.13768931](https://doi.org/10.5281/zenodo.13768931)
<br>


<br>


## License and citation

This project is licensed under the GPL 3.0-License - see the [LICENSE](LICENSE) file for details.

Please cite this repository as: Liagre, E.B.K. (2024). Medial Epicondyle Cropping Protocol. https://github.com/ElleLiagre/medial-epicondyle-cropping-protocol 
