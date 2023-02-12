

# ti_dbcti
A docker container for dbcti, mainly for cell trajectory inference analysis of single cell RNA-seq data

## dbcti repository
```shell
https://github.com/tianlt/dbcti
```

## pull container
```shell
docker pull renjun0324/ti_dbcti
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
                               
ti_dbcti = create_ti_method_container("renjun0324/ti_dbcti")
model = infer_trajectories(dataset_wrap, 
			    ti_dbcti(), 
			    parameters = list(tree_method="epg"), # or ppt
			    verbose = TRUE, return_verbose = TRUE, debug = FALSE)
model$model = map(model$model, add_cell_waypoints)
metric <- map_dfr(model$model,
                  dyneval::calculate_metrics,
                  dataset = dataset,
                  metrics = c("featureimp_wcor", "him",  "F1_branches", "correlation")) 
```                        

