

# ti_scfates_tree
A docker container for scFates, primarily for cell trajectory inference analysis of single-cell RNA-seq data.

The container was edited to follow the main scFates process, retaining some of the parameter settings, especially the tree learning methods, including 'Elastic Principal Graph (epg)' and 'simpleppt (ppt)'. Notably, we did not keep cellrank-related functions in it.


## scFates repository and help document
```shell
https://github.com/LouisFaure/scFates

https://scfates.readthedocs.io/en/latest/

```

## pull container
```shell
docker pull renjun0324/ti_scfates_tree
```

## quick start
```r

library(dynwrap)
library(dynmethods)
library(dyntoy)
library(tidyverse)
library(purrr)
library(dyno)

data("fibroblast_reprogramming_treutlein")

dataset <- wrap_expression(
  counts = fibroblast_reprogramming_treutlein$counts,
  expression = fibroblast_reprogramming_treutlein$expression
)
                               
ti_scfates_tree = create_ti_method_container("renjun0324/ti_scfates_tree:v0.1.0")
model = infer_trajectories(dataset_wrap, 
			    ti_scfates_tree(), 
			    parameters = list(tree_method="ppt"), # or epg
			    verbose = TRUE, return_verbose = TRUE, debug = FALSE)
model$model = map(model$model, add_cell_waypoints)
metric <- map_dfr(model$model,
                  dyneval::calculate_metrics,
                  dataset = dataset,
                  metrics = c("featureimp_wcor", "him",  "F1_branches", "correlation")) 
```                        

