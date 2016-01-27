# honfleuR

## Goal

Extension of `seurat` package (https://github.com/satijalab/seurat) with
following purposes:

- incorporate new schemes/algorithms to improve performance and expand
  capability.
- greatly improves running speed.
- improves codes readability by following popular [Hadley
  style](http://adv-r.had.co.nz/Style.html) to be more friendly for developers.
- workspace before I submit pull request, and hopefully it will be merged into
  `seurat` package.

## Install

```r
require(devtools)
install_github("Puriney/honfleuR", subdir = "pkg")
```

## How does `honfleuR` come?

`seurat` is the last name of the great impressionist painter Georges Seurat. The
oil painting "Evening, Honfleur" displayed at
[MoMA](http://www.moma.org/collection/works/79333?locale=en) left me great
impression which is particularly special thanks to its wooden frame. As this
package is derivative of original `seurat` package, `honfleuR` sounds good name
just like how the painting was drawn by Georges Seurat back in 1886.

## Updates

- Import PLSR-based (partial least squares regression) imputation strategy. See
  wiki [here](https://github.com/Puriney/honfleuR/wiki/Imputation-Schemes).
- Finish polishing and accelerating the main stream of inferring cell origins.
  See wiki [here](https://github.com/Puriney/honfleuR/wiki/Performance-enhancements-for-bimodal-distributions-estimation) for bimodal estimation and [here](https://github.com/Puriney/honfleuR/wiki/Performance-enhancements-for-mapping-cells-location-part) for mapping cells.
