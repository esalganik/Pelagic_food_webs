Olsen et al. (2026) Sea-ice ridge formation fuels Arctic pelagic food webs during the polar night 

MATLAB codes are provided to provide snow freeboard and ice draft from MOSAiC observations from the Airborne Laser Scanner and multibeam sonar mounted on a Remotely Operated Vehicle.
Codes are then provided to reproduce Figure 2ab in the main paper.

Figure requires selected ALS scans available at Hutter et al., 2023, doi:10.1594/PANGAEA.950896 and ROV sonar surveys available at Anhaus et al., 2024, doi:10.1594/PANGAEA.951077
Contour plot colors follow recommendations of scientifically derived color maps from Crameri et al., 2020, doi:10.1038/s41467-020-19160-7.

Repository contents:
- Fig2ab.png - reproduced figure
- Olsen_2026_Fig2ab.m - script to process raw data and reproduce Fig2ab.png
- PS122_2_19_27_20200107_ROV_MULTIBEAM_v1_new.tab - multibeam sonar dataset from Anhaus et al., 2024
- README.md - repository description
- broc.mat - colormap broc from Crameri et al., 2020
- colocated_freeboard_draft_20200107.nc - processed and exported data to reproduce the figure Fig2ab.png
