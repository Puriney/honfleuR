# honfleuR

## Goal

Refined `seurat` package (https://github.com/satijalab/seurat) with following
purposes:

- improves codes readability by following [Hadley
  style](http://adv-r.had.co.nz/Style.html).
- improves running speed by using vectorization as much as possible.
- incorporate new schemes/algorithms to improve performance.

**Hopefully it will be merged into `seurat` package**.

## Why `honfleuR`

`seurat` is the last name of the great impressionist painter Georges Seurat. The
paining "Evening, Honfleur" displayed at
[MoMA](http://www.moma.org/collection/works/79333?locale=en) left me great
impression which is particularly special thanks to the wooden frame. As this
package is derivative of original `seurat` package, `honfleuR` sounds good name
just like how the painting was drawn by Georges Seurat back in 1886.

## News

- Import PLSR-based (partial linear squared regression) imputation strategy
- Finish polishing and accelerating the main stream of inferring cell origins
