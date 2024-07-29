# Vertex-wise-LMM-analysis

This repository contains the code used for the diffusion MRI (dMRI) analysis reported in the article "Longitudinal Multi-Tensor Analysis of Neocortical Microstructure in an Animal Model of Cortical Dysplasia" [doi: 10.1101/2024.07.09.602800](https://www.biorxiv.org/content/10.1101/2024.07.09.602800v2.article-info).

## Overview

For the analysis, we developed a curvilinear coordinate system called grid_lines, which serves as a common anatomical framework for the cortex across all subjects. The grid_lines system consists of fifty lines that connect the pial surface to the white/gray matter boundary, covering the entire cortical ribbon. Each line is divided into ten equidistant points, which are used to sample the diffusion metrics. In total, we sample 500 vertices per hemisphere.

## Analysis Workflow

The analysis is carried out using a series of R scripts, following these steps:

1) Data Preparation: Tidy the data into the correct format.
2) Vertex-wise Linear Mixed-Effects Model (LMM) for each vertex.
3) Effect Size Calculation by Cohen's f for each vertex.
4) Cluster Inference Analysis.
5) Data Visualization.
