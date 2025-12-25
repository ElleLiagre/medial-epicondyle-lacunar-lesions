# Microtrauma at the Humeral Medial Epicondyle: Quantifying Lacunar Lesions

This repository contains all scripts and protocols used in the study "Microtrauma at the Humeral Medial Epicondyle: Quantifying Lacunar Lesions to Reconstruct Activity Patterns in Past Populations."

<br> 

**Author and contact details**  
Author &nbsp;&nbsp;&nbsp; *Elle B. K. Liagre*  
Email &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  *elle.liagre@u-bordeaux.fr*  
ORCID <span>&nbsp;&nbsp;&nbsp;</span>   *[https://orcid.org/0000-0002-8993-3266](https://orcid.org/0000-0002-8993-3266)*

<br>


**Table of Contents**

- [Description](#description)
- [Dependencies and libraries](#dependencies-and-libraries)
- [License and Citation](#license-and-citation)

<br>

---

## Description

This project presents a novel 3D quantitative approach for analyzing lacunar lesions at the anterior medial collateral ligament attachment site on the humeral medial epicondyle. The methodology enables objective, reproducible assessment of skeletal microtrauma in archaeological populations.

This protocol is presented in the study by **Liagre EBK, Knüsel CJ, Villotte S. Microtrauma at the Humeral Medial Epicondyle: Quantifying Lacunar Lesions to Reconstruct Activity Patterns in Past Populations (under review).** Methodological choices and further information can be found there.
<br>
<br>


**Repository Structure**
```
├── Protocol/
│   ├── Application/
│   │   ├── extraction_protocol.py                # Python script for 3D surface extraction
│   │   ├── screening_protocol_application.R      # R script for lesion screening
│   │   └── buffer_mesh_hull.RDS                  # Buffer mesh used during lesion screening
│   └── Development
│       └── screening_protocol_development.R      # R script for developing the screening protocol
├── Analysis/
│   ├── lesion_analysis.R                         # R script for statistical analysis of lesion morphology
│   ├── reliability_study.R                       # R script for intra-observer reliability analysis
│   └── volume_calculation.py                     # Python script for volume measurements of 3D meshes
├── LICENSE
└── README.md
```
<br>

## Dependencies and libraries
This project was built with the following software and library versions. 

**Software and Extensions**
- R version (>= 4.3.2)
  - Required R packages are mentioned in the individual scripts
- Python (>= 3.13.2)
  - PyMeshLab (>= 2023.12.post2)

**Citations**

- R Core Team. 2023. _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, released. https://www.R-project.org/.
- Muntoni, A., & Cignoni, P. (2024). *PyMeshLab: PyMeshLab v2023.12.post2*. Zenodo. [doi: 10.5281/zenodo.13768931](https://doi.org/10.5281/zenodo.13768931)
<br>


## License and citation

This project is licensed under the GPL 3.0-License - see the [LICENSE](LICENSE) file for details.

Please cite this repository as: Liagre, E.B.K. (2025). Medial Epicondyle Lacunar Lesions. https://github.com/ElleLiagre/medial-epicondyle-lacunar-lesions 
