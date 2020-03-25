
<!-- README.md is generated from README.Rmd. Please edit that file -->
Inference, Learning, and Optimization on Grassmann manifold
===========================================================

<!-- badges: start -->
<!-- badges: end -->
Grassmannian is a set of linear subspaces, which forms a Riemannian manifold. We provide algorithms for statistical inference, optimization, and learning over the Grassmann manifold.

Installation
------------

You can install the released version of RiemGrassmann from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RiemGrassmann")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kyoustat/RiemGrassmann")
```

Available Functions
-------------------

| function      | description                                 |
|---------------|---------------------------------------------|
| `gr.hclust`   | Hierarchical clustering.                    |
| `gr.kmedoids` | k-Medoids clustering.                       |
| `gr.mean`     | Frechet mean and variation.                 |
| `gr.pdist`    | Pairwise distance for Grassmann-valued data |
| `gr.pdist2`   | Pairwise distance between two sets of data  |
