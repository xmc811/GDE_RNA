# GDE_RNA - R Shiny Visualization App for RNA-seq Data

This repository hosts the code for the RNA-seq visualization project at The University of Texas MD Anderson Cancer Center. 

It is only for demonstration purpose and not under active development. However, it may serve as a template for RNA-seq-based visualization applications, and comments are always welcomed.

---

### To Run the Application Locally Using Docker

It is recommended to use Docker to run the R Shiny application package `rocker_gde_rna:1.0`:

```
docker pull ghcr.io/xmc811/rocker_gde_rna:1.0
docker run --rm -p 3838:3838 ghcr.io/xmc811/rocker_gde_rna:1.0
```

You should be able to access the application at `localhost:3838`.

You may also clone the repository and run it within RStudio. However, it requires installation of additional R packages including [`xmcutil`](https://github.com/xmc811/xmcutil).
