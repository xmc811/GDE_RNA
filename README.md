# GDE_RNA - R Shiny Visualization App for RNA-seq Data

This repository hosts the code for the RNA-seq visualization project at The University of Texas MD Anderson Cancer Center. 

---

## Instructions for Developers

#### 1. Overview

This application can be further developed on top of current files. It depends on several R Bionconductor packages and another R package specifically written
for bulk RNA-seq visualization [`xmcutil`](https://github.com/xmc811/xmcutil).

#### 2. Containerization

The final Docker image is built upon the `rocker_gde_base:latest` image with additional configurations in the [`Dockerfile`](https://github.com/xmc811/GDE_RNA/blob/master/Dockerfile). `rocker_gde_base:latest` is an intermediate docker image used for facilitating the CI/CD process.
During each development cycle, if changes are limited to the R Shiny App (instead of the underlying R packages), you can directly build the final image at the main directory.

The intermediate docker image `rocker_gde_base:latest` can be build upon the publicly available image [`rocker/shiny-verse:latest`](https://hub.docker.com/r/rocker/shiny-verse). The `Dockerfile` and the `.R` file for R package installations can be found [here](https://github.com/xmc811/GDE_RNA/tree/master/Rocker_GDE_base). If there are changes in these R packages, the `rocker_gde_base:latest` image needs to be re-built.

(Note: the installations of these R packages may take a long time. To increase the development efficiency, it is recommended to split the process into multiple steps, thus increasing the "layers" of the containerization.


## Instructions for Users

### To Run the Application Locally Using Docker

It is recommended to use Docker Engine to run the R Shiny application package `rocker_gde_rna:1.0`:

```
docker pull ghcr.io/xmc811/rocker_gde_rna:1.0
docker run --rm -p 3838:3838 ghcr.io/xmc811/rocker_gde_rna:1.0
```

You should be able to access the application at `localhost:3838`.

You may also clone the repository and run it within RStudio. However, it requires installation of additional R packages including [`xmcutil`](https://github.com/xmc811/xmcutil). See the Section **Instructions for Developers** for further details.

---

## User Interface

As shown in the figure below, there are four top navigation tabs. The first two tabs contain core functionalities of the app, while the last two are merely placeholders serving as templates for further development.

- Data Input & Exploration - includes data uploading, summary, sample-level data exploration, and differential gene expression (DGE) analysis.
- DGE Visualization - visualizes the DGE results, including MA plots, volcano plots, boxplots, sample-gene heatmaps, and Gene-Set Enrichment Analysis (GSEA).
- Help - a placeholder tab for embedding the `.pdf` help file.
- About - a placeholder tab for showing addtional information.

![Interface](./pics/interface.png)

---

## File Upload

The file upload section contains three options. By default, the option "Example" is selected.

- Example - sample RNA-Seq data with 10 samples.
- Upload - users can upload HTSeq-count output files, with one file per sample.
- Select - a placeholder for further development.

