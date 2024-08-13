# Vertex-wise-LMM-analysis

This repository contains the code used for the diffusion MRI (dMRI) analysis reported in the article "Longitudinal Multi-Tensor Analysis of Neocortical Microstructure in an Animal Model of Cortical Dysplasia" [doi: 10.1101/2024.07.09.602800](https://www.biorxiv.org/content/10.1101/2024.07.09.602800v2.article-info).

## Overview

For the analysis, we developed a curvilinear coordinate system called grid_lines, which provides a consistent anatomical framework for the cortex across all subjects. The grid_lines system comprises fifty lines that span from the pial surface to the white/gray matter boundary, encompassing the entire cortical ribbon. Each line is subdivided into ten equidistant points, which are used to sample diffusion metrics. Consequently, we obtain a total of 500 sampled vertices per hemisphere. 

The figure below represents the general workflow of this analysis. The upper panel illustrates the dMRI process, from raw images to adjusting the grid_lines, and how we sample our dMRI methods to each vertex. The lower panel outlines the statistical steps we applied using the R code provided in this repository.

_______

![Figura_methods](https://github.com/user-attachments/assets/d3ffb095-e00f-467e-b387-e679e218a391)

## Analysis Workflow

The analysis is carried out using a series of R scripts, following these steps:

2) Linear Mixed-Effects Model (LMM) and Effect Size by Cohen's f^2 for each vertex.
3) Permutation of LMM for each vertex.
4) Cluster inference analysis.
5) Data visualization.
