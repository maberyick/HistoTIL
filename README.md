# feature extrction of different combination of Tumor-Infiltrating Lymphocytes (TIL) and non-TIL (HistoTIL)
Morphological feature tools consisting of description of tumor-infiltrating lymphocytes (TILs) density and spatial arrangement on H&E-stained biopsy samples. Also nuclei features of quantification invloving shape and texture of non-TIL cells (tumor cells) and contextual description of similarities and differences of neighboring cells. These features are extracted across different tissue compositions (Epithelium, stroma and boundary between them).

## Input type
1. The input are .mat files, containing H&E-stained patches of size 2000x2000 generated from a Whole Slide Image (WSI).
2. A .mat file containing a binary mask identifying by a binary 1, all that are identified as nuclei (cells and membrane). A binary 0 is the background.
3. (Look here https://github.com/maberyick/LympDetect_HE_ML_V2)
4. A deep learning mask identifying with a binary 1, the epithelium tissue.
(Look here https://github.com/maberyick/EpiStroma_HE_DL_V2)

## Run
To run the script, follow the directions on the file, feat_extraction_main.m
The file will ask for,
- folder_matpatches = folder of the extracted H&E-stained .mat patches
- folder_pyepistroma = folder of the extracted epithelium masks
- folder_matcellmask = folder of the generated binary cell masks
- folder_savepath = fodler of saving the output of the tool

## Feature description
**Nuclei (n=100):** Nuclei features involve the morphological quantification of the non-TIL nuclei involving shape, color (intensity) and texture.

**Contextual (n=87):** Based on the nuclei features, different quantification metrics are calculated around each nucleus, such as how many, different or similar are, compared to the central nuclei. These metrics include, quantity of cells, their shape, area, eccentricity, number of surrounding lymphocytes.

**Density TIL(denTIL) (n=19):** It involves the density of Tumor-infiltrating Lymphocytes involving different specific metrics such as number of TILs per area, ratio of TILs and tumor cells, ratio of TIL cluster over tumor cell cluster.

**Spatial TIL (spaTIL) (n=85):** Similar to density, spatial arrangement of TILs is quantified using distinctive graph-based metrics to find niche clusters of TILs surrounding tumor clusters and their number, convex hull shape, proximity is calculated.
